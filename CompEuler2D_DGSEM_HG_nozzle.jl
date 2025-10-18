using Trixi
using OrdinaryDiffEq
using StaticArrays
using OrdinaryDiffEqLowStorageRK
using Trixi2Vtk
using LinearAlgebra
using Base.Threads
println("Julia is using $(Threads.nthreads()) threads.")

T_end=0.1 #Change if needed to get a temporally stable jet
my_poldeg=3

println("Checkpoint 1 nach T_end und my_poldeg")

#equations (Compressible Euler Eqs of argon, 5/3 is a sufficiently good approx of the ratio of specific heats gamma for a monoatomic gas)
equations = CompressibleEulerEquations2D(5/3)

println("Checkpoint 2 nach equations")

indicator = IndicatorHennemannGassner(equations, LobattoLegendreBasis(my_poldeg); alpha_max = 0.8, alpha_min = 0.001, alpha_smooth = true, variable = Trixi.density_pressure)


println("Checkpoint 3 nach Indicator")

my_VolumeIntegral= VolumeIntegralShockCapturingHG(volume_flux_dg = flux_ranocha_turbo, volume_flux_fv = flux_hllc, indicator)
#volume_flux=flux_hllc
#my_VolumeIntegral= VolumeIntegralFluxDifferencing(volume_flux)

println("Checkpoint 4 nach my_VolumeIntegral")

# Lobatto-Legendre basis for polynomials of degree my_poldeg
solver = DGSEM(polydeg=my_poldeg, surface_flux=flux_hllc, volume_integral = my_VolumeIntegral, RealT=Float64)

println("Checkpoint 5 nach solver")

#loading and working with the mesh
mesh_file = joinpath(@__DIR__, "meshes", "JetRIS_verysimp_nozzle.mesh")

println("Checkpoint 6 nach mesh_file")

mesh = UnstructuredMesh2D(mesh_file; RealT=Float64)
println("Checkpoint 7 nach mesh")


#initial condition
function my_initial_condition(x, t, equations)
    #some operating parameters
    TempJetRIS = 300.0 #K

    # x_transition_start = 0.0964 - 0.01 # Start of the transition zone
    # x_transition_end = 0.0964 + 0.01   # End of the transition zone
    x_transition_start = 0.0964
    x_transition_end = 0.0964


    pJetRIS = 8000.0 #Pa (=80 mbar)

    massnumberA_Argon = 0.039948 #kg/mol (=39.948 g/mol)
    gasconstantR = 8.3145 #J/(mol*K)
    rhoJetRIS = massnumberA_Argon * pJetRIS / (gasconstantR * TempJetRIS)

    v1_IC  = 0.0
    v2_IC  = 0.0

	# Calculate conserved variables
	rho_E_IC = pJetRIS / (equations.gamma - 1.0) + 0.5 * rhoJetRIS * (v1_IC^2.0 + v2_IC^2.0) # is p / (gamma - 1) because 0.5 * rho * (v1^2 + v2^2) is zero
	return SVector(rhoJetRIS, rhoJetRIS * v1_IC, rhoJetRIS * v2_IC, rho_E_IC)
end

println("Checkpoint 8 nach my_initial_condition")

#Boundary conditions
#inflow
@inline function inflow_variables_function(x, t, equations)
   #JetRIS properties etc. for calculating inflow parameters

    # general constants
    gamma = equations.gamma
    massnumberA_Argon = 0.039948 #kg/mol (=39.948 g/mol)
    gasconstantR = 8.3145 #J/(mol*K)
    cp_Argon = gamma * gasconstantR / ((gamma - 1) * massnumberA_Argon)

    #stagnation temp
    Temp_stagn = 300.0 #K

    #JetRIS measurements, claculation of volumetric flow rate etc. and calculation of x velocity
    diamINLETJetRIS = 0.0006 #m (=0.6 mm)
    volumetricflowrateQ =52*diamINLETJetRIS^2.0*sqrt(Temp_stagn/massnumberA_Argon) #Formula from Yu. Kudryavtsev et al 2012 in SI units
    diamCELLJetRIS = 0.097 #m (= 97 mm)
    CrosssectionalareaJetRIS = pi*(diamCELLJetRIS/2.0)^2.0

 
    v1_ramp_up  = volumetricflowrateQ / CrosssectionalareaJetRIS
    v2_ramp_up  = 0.0

    #calculation of inflow temp
    Temp_inflow = Temp_stagn-v1_ramp_up^2/2cp_Argon

    #ratio of spec heats and calc of speed of sound and mach number
    
    speedofsound = sqrt(gamma*gasconstantR*Temp_inflow/massnumberA_Argon)
    Mach = v1_ramp_up/speedofsound
    p_stagn   = 8000.0 #Pa (=80 mbar)

    #stagnation pressure and calculation of inflow pressure
    p_inflow = p_stagn/((1+((gamma-1)/2)*Mach^2)^(gamma/(gamma-1)))

    #calculation of inflow density
    rho_inflow = massnumberA_Argon * p_inflow / (gasconstantR * Temp_inflow)
    ramp_up_time = 1.0e-2

     if t < ramp_up_time
        # ramp function
        v1_current = v1_ramp_up * (0.5 * (1.0 + tanh(3.0 * (t - ramp_up_time / 2.0) / ramp_up_time)))
        v2_current = v2_ramp_up
    else
        v1_current = v1_ramp_up
        v2_current = v2_ramp_up
    end

    gamma = equations.gamma
    rho_E_inflow = p_inflow / (gamma - 1.0) + 0.5 * rho_inflow * (v1_current^2.0 + v2_current^2.0)
    return SVector(rho_inflow, rho_inflow * v1_current, rho_inflow * v2_current, rho_E_inflow)
end

@inline function inflow_boundary_condition(u_inner, normal_direction::AbstractVector, x, t, surface_flux_function, equations::CompressibleEulerEquations2D)
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

        p_local = pressure(inflow_variables_function(x, t, equations), equations)
    end

    # Create the `u_surface` solution state where the local pressure is possibly set from an external value
    prim = SVector(rho_local, v_normal, v_tangent, p_local)
    u_boundary = prim2cons(prim, equations)
    u_surface = Trixi.rotate_from_x(u_boundary, normal, equations)

    # Compute the flux using the appropriate mixture of internal / external solution states
    return flux(u_surface, normal_direction, equations)
end


println("Checkpoint 9 nach inflow_state_function")


#outflow
@inline function outflow_variables_function(x, t, equations)
    gamma = equations.gamma
    massnumberA_Argon = 39.948e-3 #kg/mol
    gasconstantR = 8.3145 #J/(mol*K)

    cp_Argon = gamma * gasconstantR / ((gamma - 1) * massnumberA_Argon)
    Temp_stagn = 300.0 
    
    #p_stagn    = 18 #Pa
    #p_stagn    = 16 #Pa 
    #p_stagn    = 20 #Pa 
    #p_stagn    = 14
    #p_stagn    = 1 #Pa 
    #p_stagn    = 4 #Pa
    p_stagn    = 5 #Pa 
    #p_stagn    = 6 #Pa 
    #p_stagn    = 7 #Pa 
    #p_stagn    = 10 #Pa 
    #Uncomment if needed

    # stagnation pressure and calculation of outflow pressure
    massnumberA_Argon = 39.948e-3 # kg/mol
    gasconstantR = 8.3145 # J/(mol*K)
    

    pumpvol_outflow = 1.3 # m^3/s; max 1.3 m^3/s
    outflow_diam = 0.197 # m 
    outflow_surface = pi*(outflow_diam/2)^2

    v1_outflow = pumpvol_outflow/outflow_surface #(https://www.sigmaaldrich.com/DE/de/technical-documents/technical-article/protein-biology/protein-purification/converting-flow-velocity-volumetric-flow-rates?srsltid=AfmBOoq6nttgSgCeh4Y4ba-prXuYncg3by7TAPL11OQw2zL0eupqyARD)
    v2_outflow = 0.0

    Temp_outflow = Temp_stagn-v1_outflow^2/2cp_Argon #Temp_stagn-(sqrt(v1^2+v2^2))^2/2cp_argon but v2 is zero

    speedofsound = sqrt(gamma*gasconstantR*Temp_outflow/massnumberA_Argon)
    Mach = v1_outflow/speedofsound # sqrt(v1^2+v2^2)/speedofsound, but v2 is zero
    p_outflow = p_stagn/((1+((gamma-1)/2)*Mach^2)^(gamma/(gamma-1)))

    rho_outflow = p_outflow * massnumberA_Argon / (gasconstantR * Temp_outflow)

    gamma = equations.gamma
    rho_E_outflow = p_outflow / (gamma - 1.0) + 0.5 * rho_outflow * (v1_outflow^2.0 + v2_outflow^2.0)
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
    :l_vwall_inflow_symmetric_R      => inflow_boundary_condition,  # inflow
    :r_vwall_outflow_symmetric_R     => outflow_boundary_condition, # outflow
    :r_u_hwall_upper_R               => boundary_condition_slip_wall,
    :r_u_vwall_rightorange_R         => boundary_condition_slip_wall,
    :r_u_vwall_nozzlebackside_first_R => boundary_condition_slip_wall,
    :nozzle_spline_R                 => boundary_condition_slip_wall,
    :r_u_diagwall_funnel_first_R     => boundary_condition_slip_wall,
    :r_u_circwall_funnel_second_R    => boundary_condition_slip_wall,
    :r_u_hwall_funnel_third_R        => boundary_condition_slip_wall,
    :r_u_diagwall_funnel_fourth_R    => boundary_condition_slip_wall,
    :r_u_circwall_funnel_fifth_R     => boundary_condition_slip_wall,
    :r_u_vwall_funnel_sixth_R        => boundary_condition_slip_wall,
    :r_u_hwall_funnel_seventh_R      => boundary_condition_slip_wall,
    :r_u_vwall_leftorange_R          => boundary_condition_slip_wall,
    :l_u_hwall_upper_R               => boundary_condition_slip_wall,
    :l_vwall_inflow_symmetric        => inflow_boundary_condition,  # inflow
    :r_vwall_outflow_symmetric       => outflow_boundary_condition, # outflow
    :r_u_hwall_upper                 => boundary_condition_slip_wall,
    :r_u_vwall_rightorange           => boundary_condition_slip_wall,
    :r_u_vwall_nozzlebackside_first  => boundary_condition_slip_wall,
    :nozzle_spline                   => boundary_condition_slip_wall,
    :r_u_diagwall_funnel_first       => boundary_condition_slip_wall,
    :r_u_circwall_funnel_second      => boundary_condition_slip_wall,
    :r_u_hwall_funnel_third          => boundary_condition_slip_wall,
    :r_u_diagwall_funnel_fourth      => boundary_condition_slip_wall,
    :r_u_circwall_funnel_fifth       => boundary_condition_slip_wall,
    :r_u_vwall_funnel_sixth          => boundary_condition_slip_wall,
    :r_u_hwall_funnel_seventh        => boundary_condition_slip_wall,
    :r_u_vwall_leftorange            => boundary_condition_slip_wall,
    :l_u_hwall_upper                 => boundary_condition_slip_wall
)



println("Checkpoint 10 nach BCs")


#callback saves H5 files
solution_callback = SaveSolutionCallback(dt = 0.001, save_initial_solution = true, save_final_solution = true, solution_variables = cons2prim, output_directory=joinpath(@__DIR__, "h5_raw_data"))
alive_callback = AliveCallback(analysis_interval = 0, alive_interval = 100)

callbacks = CallbackSet(solution_callback, alive_callback)

println("Checkpoint 11 nach callback")

#Bundeling, method of lines, solving the ODEs
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

