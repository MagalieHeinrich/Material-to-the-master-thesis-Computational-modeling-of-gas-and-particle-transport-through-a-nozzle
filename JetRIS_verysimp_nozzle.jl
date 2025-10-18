using GLMakie, HOHQMesh

JetRIS_verysimp_nozzle = newProject("JetRIS_verysimp_nozzle", joinpath(@__DIR__, "meshes"))
setPolynomialOrder!(JetRIS_verysimp_nozzle, 1)
setMeshFileFormat!(JetRIS_verysimp_nozzle, "ISM-V2")

scaling_factor = 1.0 / 1000.0 



lower_left = [-0.221, -0.005, 0.0] 
#Comment the 18 mm option and uncomment an other option for a finer mesh
spacing = [0.018, 0.018, 0.0]
num_intervals = [70, 20, 0] 
#spacing = [0.009, 0.009, 0.0]
#num_intervals = [140, 40, 0] 
#spacing = [0.0045, 0.0045, 0.0]
#num_intervals = [280, 80, 0] 
#spacing = [0.00225, 0.00225, 0.0]
#num_intervals = [560, 160, 0] 

addBackgroundGrid!(JetRIS_verysimp_nozzle, lower_left, spacing, num_intervals)
plotProject!(JetRIS_verysimp_nozzle, GRID)

b=1*0.108# Change this for a wider y-diameter
a=0.399# Change this for a longer right chamber
l_u_hwall_upper = newEndPointsLineCurve("l_u_hwall_upper",[0.085,0.126+b,0.000],[-0.215,0.126+b,0.000])

l_vwall_inflow_symmetric = newEndPointsLineCurve("l_vwall_inflow_symmetric",[-0.215,0.126+b,0.000],[-0.215,0.000,0.000])

sym_inlet_to_funnel_start = newEndPointsLineCurve(":symmetry", [-0.215, 0.0, 0.0], [0.085, 0.0, 0.0])
sym_funnel_to_nozzle_throat = newEndPointsLineCurve(":symmetry", [0.085, 0.0, 0.0], [0.0934, 0.0, 0.0])
sym_nozzle_mid = newEndPointsLineCurve(":symmetry", [0.0934, 0.0, 0.0], [0.1144, 0.0, 0.0])
sym_nozzle_exit_to_outflow = newEndPointsLineCurve(":symmetry", [0.1144, 0.0, 0.0], [0.415+a, 0.0, 0.0])

r_vwall_outflow_symmetric = newEndPointsLineCurve("r_vwall_outflow_symmetric",[0.415+a,0.000,0.000],[0.415+a,0.126+b,0.000])
r_u_hwall_upper = newEndPointsLineCurve("r_u_hwall_upper",[0.415+a,0.126+b,0.000],[0.105,0.126+b,0.000])
r_u_vwall_rightorange = newEndPointsLineCurve("r_u_vwall_rightorange",[0.105,0.126+b,0.000],[0.105,0.0113,0.000])
r_u_vwall_nozzlebackside_first = newEndPointsLineCurve("r_u_vwall_nozzlebackside_first",[0.105,0.0113,0.000],[0.1144,0.0113,0.000])
r_u_hwall_nozzlebackside_second = newEndPointsLineCurve("r_u_hwall_nozzlebackside_second",[0.1144,0.0113,0.000],[0.1144,0.0093,0.000])

nozzle_spline = newSplineCurve("nozzle_spline", joinpath(@__DIR__, "meshes", "splinedata_nozzle.txt"))#spline from spline data

r_u_diagwall_funnel_first = newEndPointsLineCurve("r_u_diagwall_funnel_first",[0.0934,0.003,0.000],[0.088229+0.005*sind(90-44.749),0.001+0.005*cosd(90-44.749),0.000])
r_u_circwall_funnel_second = newCircularArcCurve("r_u_circwall_funnel_second",[0.088229,0.001,0.000],0.005,44.749,90.000,"degrees")
r_u_hwall_funnel_third = newEndPointsLineCurve("r_u_hwall_funnel_third",[0.088229,0.006,0.000],[0.086,0.006,0.000])
r_u_diagwall_funnel_fourth = newEndPointsLineCurve("r_u_diagwall_funnel_fourth",[0.086,0.006,0.000],[0.071-cosd(247.339-180)*0.001,0.013346-sind(247.339-180)*0.001,0.000])
r_u_circwall_funnel_fifth = newCircularArcCurve("r_u_circwall_funnel_fifth",[0.071,0.013346,0.000],0.001,247.339,180.000,"degrees")
r_u_vwall_funnel_sixth = newEndPointsLineCurve("r_u_vwall_funnel_sixth",[0.070,0.013346,0.000],[0.070,0.0229,0.000])
r_u_hwall_funnel_seventh = newEndPointsLineCurve("r_u_hwall_funnel_seventh",[0.070,0.0229,0.000],[0.085,0.0229,0.000])
r_u_vwall_leftorange = newEndPointsLineCurve("r_u_vwall_leftorange",[0.085,0.0229,0.000],[0.085,0.126+b,0.000])

#adding the curves
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,l_u_hwall_upper)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,l_vwall_inflow_symmetric)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,sym_inlet_to_funnel_start)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,sym_funnel_to_nozzle_throat)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,sym_nozzle_exit_to_outflow)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,sym_nozzle_mid)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_vwall_outflow_symmetric)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_hwall_upper)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_vwall_rightorange)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_vwall_nozzlebackside_first)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_hwall_nozzlebackside_second)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,nozzle_spline)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_diagwall_funnel_first)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_circwall_funnel_second)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_hwall_funnel_third)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_diagwall_funnel_fourth)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_circwall_funnel_fifth)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_vwall_funnel_sixth)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_hwall_funnel_seventh)
addCurveToOuterBoundary!(JetRIS_verysimp_nozzle,r_u_vwall_leftorange)

#refinement
h=0.500
w=0.200

ref_u_funnelleftside=newRefinementLine("line", "smooth",  [ 70.000*scaling_factor,22.900*scaling_factor, 0.000],[ 70.000*scaling_factor,13.350*scaling_factor, 0.000],h*scaling_factor, w*scaling_factor)
ref_middleline=newRefinementLine("line","smooth",[0.095,0.0,0.0],[0.098,0.0,0.0],h*scaling_factor, w*scaling_factor)
ref_r_u_diagwall_nozzle_first = newRefinementLine("line", "smooth",[0.1144,0.0093,0.000],[0.107060,0.005540,0.000],h*scaling_factor, w*scaling_factor)
ref_r_u_diagwall_nozzle_second = newRefinementLine("line", "smooth",[0.107060,0.005540,0.000],[0.102525,0.003217,0.000],h*scaling_factor, w*scaling_factor)
ref_r_u_diagwall_nozzle_third = newRefinementLine("line", "smooth",[0.102525,0.003217,0.000],[0.098050,0.000924,0.000],h*scaling_factor, w*scaling_factor)
ref_u_noz_extra=newRefinementLine("line","smooth",[0.096425,0.000501,0.000],[0.0934,0.003,0.000],h*scaling_factor, w*scaling_factor)
ref_missingcurve=newRefinementLine("line","smooth",[0.098050, 0.000924, 0.000],[0.096425, 0.000501, 0.000],h*scaling_factor, w*scaling_factor)

addRefinementRegion!(JetRIS_verysimp_nozzle,ref_u_funnelleftside)
addRefinementRegion!(JetRIS_verysimp_nozzle,ref_middleline) 
addRefinementRegion!(JetRIS_verysimp_nozzle,ref_r_u_diagwall_nozzle_second)
addRefinementRegion!(JetRIS_verysimp_nozzle,ref_r_u_diagwall_nozzle_third)
addRefinementRegion!(JetRIS_verysimp_nozzle,ref_u_noz_extra)
addRefinementRegion!(JetRIS_verysimp_nozzle,ref_missingcurve)

generate_mesh(JetRIS_verysimp_nozzle)
