# CliMA Soil Modeling Code
One can use the code here, at present, to simulate Richard's equation, which describes how water flows in soil.

To run, from this directory, enter
```
> julia --project=../../../ soil_water_test.jl
```

This currently saves the final output state into a file, and creates a plot of the moisture content vs. depth. It's pretty straightfoward to save the output at a different step, all the output, or to plot the other variables.

We currently support different formulations of the matric potential and the hydraulic conductivity. To run with a different hydraulic model, replace the WF definition in ```soil_water_test.jl```.

Note that the "Havercamp formulation" for the hydraulic conductivity is paired with the van Genuchten matric potential. The default parameters are for Yolo light clay. This corresponds to the example in Bonan's book, Ch 8, Figure 8.9. That is, the Havercamp formulation for Yolo light clay is specified by:
```
WF = waterfunctions(
    hydraulic_cond = Havercamp{FT}(),
    matric_pot = vanGenuchten{FT}()
)
```

To change the defaults, or to change to the full van Genuchten formulation, do something like:
```
WF = waterfunctions(
    hydraulic_cond = vanGenuchten{FT}(n = 2.08 ,α = 2.0),
    matric_pot = vanGenuchten{FT}(n = 2.08 ,α = 2.0)
)
```
There is also the option to use the Brooks and Corey hydraulic model, but it isn't likely that we will use this in CliMA.

There are currently issues with numerical stability of solutions that are being addressed.
