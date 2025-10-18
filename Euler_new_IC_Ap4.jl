using Trixi
using OrdinaryDiffEq
using StaticArrays
using OrdinaryDiffEqLowStorageRK
using Trixi2Vtk
using LinearAlgebra
using Base.Threads
println("Julia is using $(Threads.nthreads()) threads.")

T_end=0.2 #Change if needed to get a temporally stable jet
my_poldeg=3

println("Checkpoint 1 after T_end und my_poldeg")

#equations (Compressible Euler Eqs of argon, 5/3 is a sufficiently good approx of the ratio of specific heats gamma for a monoatomic gas)
equations = CompressibleEulerEquations2D(5/3)

println("Checkpoint 2 after equations")

indicator = IndicatorHennemannGassner(equations, LobattoLegendreBasis(my_poldeg); alpha_max = 1.0, alpha_min = 0.001, alpha_smooth = true, variable = Trixi.density_pressure)

println("Checkpoint 3 after Indicator")

my_VolumeIntegral= VolumeIntegralShockCapturingHG(volume_flux_dg = flux_ranocha_turbo, volume_flux_fv = flux_hllc, indicator)
#my_VolumeIntegral= VolumeIntegralShockCapturingHG(volume_flux_dg = flux_ranocha_turbo, volume_flux_fv = flux_lax_friedrichs, indicator) #Uncomment if more dissipation is needed

#volume_flux = flux_ranocha_turbo
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux) #Uncomment if shock capturing is unwanted


println("Checkpoint 4 after my_VolumeIntegral")

# Lobatto-Legendre basis for polynomials of degree my_poldeg
solver = DGSEM(polydeg=my_poldeg, surface_flux=flux_hllc, volume_integral = my_VolumeIntegral, RealT=Float64)
#solver = DGSEM(polydeg=my_poldeg, surface_flux=flux_lax_friedrichs, volume_integral = my_VolumeIntegral, RealT=Float64) #Uncomment if more dissipation is needed


println("Checkpoint 5 after solver")

#loading and working with the mesh
mesh_file = joinpath(@__DIR__, "meshes", "JetRIS_verysimp_ap5.mesh")

println("Checkpoint 6 after mesh_file")

mesh = UnstructuredMesh2D(mesh_file; RealT=Float64)
println("Checkpoint 7 after mesh")


#initial condition
function my_initial_condition(x, t, equations)
    
    # some constants
    gamma = equations.gamma
    massnumberM_m_Argon = 0.039948 #kg/mol (=39.948 g/mol)  #kg/mol https://www.ciaaw.org/argon.htm
    gasconstantR = 8.3145 #J/(mol*K) Atkins, title page
    
    #x_transition_start = 0.09 - 0.01 # Start of the transition zone
    #x_transition_end = 0.09 + 0.01   # End of the transition zone
    x_transition_start = 0.0964
    x_transition_end = 0.0964

    #stagnation pressure
    p_stagn = 8000 #Pa
    #ambient pressure
    p_amb = 5 #Pa
    
    #Mach number
    Mach = sqrt((2*((p_stagn/p_amb)^((gamma-1)/gamma))-2)/(gamma-1))

    #stagnation temperature
    Temp_stagn = 300.0 #K
    #ambient temperature
    Temp_amb = Temp_stagn/(1+((gamma-1)/2)*Mach^2)

    #speed of sound
    speedofsound = sqrt(gamma*gasconstantR*Temp_stagn/massnumberM_m_Argon)

    #stagnation density
    rho_stagn = massnumberM_m_Argon * p_stagn / (gasconstantR * Temp_stagn)
    #ambient density
    rho_amb = rho_stagn*((1+((gamma-1)/2)*Mach^2)^(-1/(gamma-1)))

    #stagnation velocities
    v1_stagn  = 0.0
    v2_stagn  = 0.0
    # ambient velocities
    v1_amb = speedofsound*Mach
    v2_amb = 0.0

	# stagnation total energy 
	rho_E_stagn = p_stagn / (gamma - 1.0) + 0.5 * rho_stagn * (v1_stagn^2 + v2_stagn^2) # is p / (gamma - 1) because 0.5 * rho * (v1^2 + v2^2) is zero
	rho_E_amb = p_amb / (gamma - 1.0) + 0.5 * rho_amb * (v1_amb^2 + v2_amb^2) # is p / (gamma - 1) because 0.5 * rho * (v1^2 + v2^2) is zero
	
    #Zones
    if x[1] < x_transition_start
        p_initial = p_stagn
        rho_initial = rho_stagn
        v1_initial = v1_stagn
        v2_initial = v2_stagn
        rho_E_initial = rho_E_stagn

    elseif x[1] > x_transition_end
        p_initial = p_amb
        rho_initial = rho_amb
        v1_initial = v1_amb
        v2_initial = v2_amb
        rho_E_initial = rho_E_amb
    else
        p_slope = (p_amb - p_stagn) / (x_transition_end - x_transition_start)
        p_initial = p_stagn + p_slope * (x[1] - x_transition_start)
        p_slope = (p_amb - p_stagn) / (x_transition_end - x_transition_start)
        p_initial = p_stagn + p_slope * (x[1] - x_transition_start)
        rho_slope = (rho_amb - rho_stagn) / (x_transition_end - x_transition_start)
        rho_initial = rho_stagn + rho_slope * (x[1] - x_transition_start)
        v1_slope = (v1_amb - v1_stagn) / (x_transition_end - x_transition_start)
        v1_initial = v1_stagn + v1_slope * (x[1] - x_transition_start)
        v2_slope = (v2_amb - v2_stagn) / (x_transition_end - x_transition_start)
        v2_initial = v2_stagn + v2_slope * (x[1] - x_transition_start)
        rho_E_slope = (rho_E_amb - rho_E_stagn) / (x_transition_end - x_transition_start)
        rho_E_initial = rho_E_stagn + rho_E_slope * (x[1] - x_transition_start)
    end
    
    return SVector(rho_initial, rho_initial * v1_initial, rho_initial * v2_initial, rho_E_initial)
end

println("Checkpoint 8 after my_initial_condition")

#Boundary conditions
#inflow
@inline function inflow_variables_function(x, t, equations)

    # some constants
    gamma = equations.gamma
    massnumberM_m_Argon = 0.039948 #kg/mol (=39.948 g/mol)  #kg/mol https://www.ciaaw.org/argon.htm
    gasconstantR = 8.3145 #J/(mol*K) Atkins, title page
    
    v2_inflow  = 0.0
    v1_inflow  = 0.0

    #inflow temp
    Temp_inflow = 300 #K

    #inflow pressure
    p_inflow = 8000 #Pa
    
    #inflow density
    rho_inflow = massnumberM_m_Argon * p_inflow / (gasconstantR * Temp_inflow) #Formeslammlung p.90,95,96, derived

    #inflow total energy 
    rho_E_inflow = p_inflow / (gamma - 1.0) + 0.5 * rho_inflow * (v1_inflow^2.0 + v2_inflow^2.0)  #ShockCap Paper S.10
    return SVector(rho_inflow, rho_inflow * v1_inflow, rho_inflow * v2_inflow, rho_E_inflow)
end

@inline function inflow_boundary_condition(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function, equations::CompressibleEulerEquations2D)
    u_boundary = inflow_variables_function(x, t, equations)

    #return surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    return flux_hllc(u_inner, u_boundary, normal_direction, equations)
end


println("Checkpoint 9 after inflow_state_function")


#outflow
@inline function outflow_variables_function(x, t, equations)
    #some constants
    gamma = equations.gamma
    massnumberM_m_Argon = 39.948e-3 #kg/mol https://www.ciaaw.org/argon.htm
    gasconstantR = 8.3145 #J/(mol*K) Atkins, title page

    #heat capacity with constant pressure per mass
    cp_Argon = gamma * gasconstantR / ((gamma - 1) * massnumberM_m_Argon)#capacité thermique massique [J/(kg K)] p.29 Bureau international des poids et mesures, SI
    
    #ambient temperature
    Temp_outflow = 300.0 #JetRIS Paper
    
    #ambient pressure
    p_outflow    = 5 #Pa JetRISPaper: 1-10 Pa

    p_amb = 5 #Pa
    p_stagn =8000 #Pa

    #computation of the outflow pressure
    speedofsound = sqrt(gamma*gasconstantR*Temp_outflow/massnumberM_m_Argon) #Atkins S. 81
    Mach = sqrt((2*((p_stagn/p_amb)^((gamma-1)/gamma))-2)/(gamma-1))

    #v1_outflow = pumpvol_outflow/outflow_surface #(https://www.sigmaaldrich.com/DE/de/technical-documents/technical-article/protein-biology/protein-purification/converting-flow-velocity-volumetric-flow-rates?srsltid=AfmBOoq6nttgSgCeh4Y4ba-prXuYncg3by7TAPL11OQw2zL0eupqyARD)
    v1_outflow = speedofsound*Mach
    v2_outflow = 0.0

    #outflow density
    rho_outflow = p_outflow * massnumberM_m_Argon / (gasconstantR * Temp_outflow) #Formeslammlung S.90,95,96, derived

    #outflow total energy 
    rho_E_outflow = p_outflow / (gamma - 1.0) + 0.5 * rho_outflow * (v1_outflow^2.0 + v2_outflow^2.0) #ShockCap Paper S.10
    return SVector(rho_outflow, rho_outflow * v1_outflow, rho_outflow * v2_outflow, rho_E_outflow)
end

@inline function outflow_boundary_condition(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function, equations::CompressibleEulerEquations2D)
    # This would be for the general case where we need to check the magnitude of the local Mach number
    norm_ = norm(normal_direction)
    # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
    normal = normal_direction / norm_

    # Rotate the internal solution state
    u_local = Trixi.rotate_to_x(u_inner, normal, equations)

    # Compute the primitive variables
    rho_local, v_normal, v_tangent, p_local = cons2prim(u_local, equations)

    # Compute local Mach number
    a_local = sqrt(equations.gamma * p_local / rho_local)
    Mach_local = abs(v_normal / a_local)
    if Mach_local <= 1.0 
 
        p_local = pressure(outflow_variables_function(x, t, equations), equations)
    end

    # Create the `u_surface` solution state where the local pressure is possibly set from an external value
    prim = SVector(rho_local, v_normal, v_tangent, p_local)
    u_boundary = prim2cons(prim, equations)
    u_surface = Trixi.rotate_from_x(u_boundary, normal, equations)

    # Compute the flux using the appropriate mixture of internal / external solution states
    return flux(u_surface, normal_direction, equations)
end


my_boundary_conditions = Dict(
    
    :l_ulw_vwall_inflow_R           => inflow_boundary_condition,  # inflow
    :l_lw_hwall_lower_R             => boundary_condition_slip_wall,
    :r_lw_vwall_leftorange_R        => boundary_condition_slip_wall,
    :r_lw_circwall_funnel_third_R   => boundary_condition_slip_wall,
    :r_lw_hwall_funnel_sixth_R      => boundary_condition_slip_wall,
    :r_lw_circwall_funnel_seventh_R => boundary_condition_slip_wall,
    :r_lw_diagwall_funnel_eighth_R  => boundary_condition_slip_wall,
    :r_lw_vwall_funnel_nineth_R     => boundary_condition_slip_wall,
    :r_lw_hwall_lower_R             => boundary_condition_slip_wall,
    :r_lwu_vwall_outflow_R          => outflow_boundary_condition, # outflow
    :l_ulw_vwall_inflow             => inflow_boundary_condition,  # inflow
    :l_lw_hwall_lower               => boundary_condition_slip_wall,
    :r_lw_vwall_leftorange          => boundary_condition_slip_wall,
    :r_lw_circwall_funnel_third     => boundary_condition_slip_wall,
    :r_lw_hwall_funnel_sixth        => boundary_condition_slip_wall,
    :r_lw_circwall_funnel_seventh   => boundary_condition_slip_wall,
    :r_lw_diagwall_funnel_eighth    => boundary_condition_slip_wall,
    :r_lw_vwall_funnel_nineth       => boundary_condition_slip_wall,
    :r_lw_hwall_lower               => boundary_condition_slip_wall,
    :r_lwu_vwall_outflow            => outflow_boundary_condition # outflow
    )


println("Checkpoint 10 after BCs")


#callbacks
#dt = 0.001
solution_callback = SaveSolutionCallback(dt = 0.0001, save_initial_solution = true, save_final_solution = true, solution_variables = cons2prim, output_directory=joinpath(@__DIR__, "h5_raw_data"))
alive_callback = AliveCallback(analysis_interval = 0, alive_interval = 100)

callbacks = CallbackSet(solution_callback, alive_callback)

println("Checkpoint 11 nach callback")

#Bundeling, method of lines, solving the ODEs
#setprecision(BigFloat, 256)
semi = SemidiscretizationHyperbolic(mesh, equations, my_initial_condition, solver; boundary_conditions = my_boundary_conditions)

println("Checkpoint 12 nach semi")

ode = semidiscretize(semi, (0.0, T_end))

println("Checkpoint 13 nach ode")

stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds = (1.0e-12, 1.0e-12),
                                                     variables = (pressure, Trixi.density))
sol = solve(ode, SSPRK43(stage_limiter!); abstol = 1.0e-5, reltol = 1.0e-5, ode_default_options()..., callback = callbacks, maxiters = 1e9)

println("Checkpoint 14 nach solve")

trixi2vtk(joinpath(@__DIR__, "h5_raw_data/solution_*.h5"), output_directory=joinpath(@__DIR__, "vtk_end_data"))

println("Checkpoint 15 nach trixi2vtk")

