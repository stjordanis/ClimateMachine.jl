# # This is a tutorial showing
# # - how to run the Soil Moisture model (Richard's Equation)
# # - how to change between a van Genuchten hydraulic model and a Havercamp hydraulic model.
# The files used are: soil_water_test.jl (the only one you should need to modify), soil_water_model_test.jl (where the boundary conditions are specified, where the balance law is defined), and where the hydraulic model functions and auxiliary functions are defined (Water/water_functions.jl, Water/hydraulic_conductivity_composable.jl, and Water/matric_potential_composable.jl).
#

# From the ClimateMachine.jl/tutorials/Land/Soil directory:
# specify the project location - should be the directory ClimateMachine.jl/

katherinedeck@KMD Soil % julia --project=../../../ soil_water_test.jl
1) Import/Export Needed Functions
2) Set up parameters of system, boundary conditions, and initial conditions...
3) Define numerical configuration for simulation...
┌ Info: Model composition
│     param_set = EarthParameterSet()
│     WF = waterfunctions{Havercamp{Float64},vanGenuchten{Float64}}(Havercamp{Float64}(1.77, 0.0359350325289575, 0.0359350325289575), vanGenuchten{Float64}(1.43, 0.3006993006993006, 2.6))
│     initialκ = #45
│     initialν = #46
│     surfaceν = #47
└     initialh = #48
┌ Info: Establishing single stack configuration for SoilMoistureModel
│     precision        = Float64
│     polynomial order = 5
│     domain_min       = 0.00 m x0.00 m x-1.00 m
│     domain_max       = 1.00 m x1.00 m x0.00 m
│     MPI ranks        = 1
│     min(Δ_horz)      = 0.12 m
└     min(Δ_vert)      = 0.01 m
[ Info: Initializing SoilMoistureModel
4) Run the simulation, and save the output...
┌ Info: Starting SoilMoistureModel
│     dt              = 6.00000e+00
│     timeend         = 43200.00
│     number of steps = 7200
└     norm(Q)         = 2.3999999999999980e-01
┌ Info: Update
│     simtime = 13740.00 / 43200.00
│     runtime = 00:00:59
└     norm(Q) = 2.5556600203997487e-01
┌ Info: Update
│     simtime = 31278.00 / 43200.00
│     runtime = 00:01:59
└     norm(Q) = 2.6681851199624107e-01
┌ Info: Finished
│     norm(Q)            = 2.7234449602564659e-01
│     norm(Q) / norm(Q₀) = 1.1347687334401950e+00
└     norm(Q) - norm(Q₀) = 3.2344496025646791e-02
keys(all_data[0]) = Any["t", "ν", "κ", "z", "h"]
#----------------------------------------------------------------------------

# This currently saves the final output state into a file, and creates a plot of the moisture content vs. depth. It's pretty straightfoward to save the output at a different step, all the output, or to plot the other variables.


# To run with a different hydraulic model, replace the WF definition in soil_water_test.jl.

# Note that the "Havercamp formulation" for the hydraulic conductivity is paired with the van Genuchten matric potential. The default parameters are for Yolo light clay. This corresponds to the example in Bonan's book, Ch 8, Figure 8.9. That is, the Havercamp formulation for Yolo light clay is specified by:
WF = waterfunctions(
    hydraulic_cond = Havercamp{FT}(),
    matric_pot = vanGenuchten{FT}()
)
#----------------------------------------------------------------------------
# To change the defaults, or to change to the full van Genuchten formulation, do something like:
WF = waterfunctions(
    hydraulic_cond = vanGenuchten{FT}(n = 2.08 ,α = 2.0),
    matric_pot = vanGenuchten{FT}(n = 2.08 ,α = 2.0)
)
#----------------------------------------------------------------------------
# You may want to use the plotting tutorial (to be made) to see how these functions depend on these parameters. The parameters above were chosen such that the profiles look somewhat similar to the Havercamp formulation for hydraulic conductivity and matric potential.  

# You also made need a smaller timestep - if you get domain errors in the numerical integrations, try a smaller dt.

# There is also the option to use the Brooks and Corey hydraulic model, but it isn't likely that we will use this in CliMA.

