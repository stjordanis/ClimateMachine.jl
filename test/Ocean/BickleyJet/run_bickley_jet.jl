#!/usr/bin/env julia --project
using Test
using ClimateMachine
ClimateMachine.init()
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.VariableTemplates
using ClimateMachine.Mesh.Grids: polynomialorders
using ClimateMachine.Ocean.BickleyJet
using ClimateMachine.Ocean.OceanProblems

using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.BalanceLaws: vars_state, Prognostic, Auxiliary
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.VTK

using MPI
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates

function run_bickley_jet(dt, nout)
    mpicomm = MPI.COMM_WORLD
    ArrayType = ClimateMachine.array_type()

    brickrange = (xrange, yrange)
    topl = BrickTopology(
        mpicomm,
        brickrange,
        periodicity = (true, true),
        boundary = ((0, 0), (0, 0)),
    )
    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = N,
    )

    problem = SimpleBox{FT}(Lˣ, Lʸ, H)

    model = BickleyJetModel{FT}(
        problem,
        NonLinearAdvectionTerm(),
        BickleyJet.ConstantViscosity{FT}(
            ν = 5e3,   # m²/s
            κ = 1e3,   # m²/s
        ),
        nothing,
        nothing;
        g = 10, # m/s²
        c = 2, # m/s
    )

    dg = DGModel(
        model,
        grid,
        RusanovNumericalFlux(),
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient(),
    )

    Q = init_ode_state(dg, FT(0); init_on_cpu = true)

    lsrk = LSRK54CarpenterKennedy(dg, Q, dt = dt, t0 = 0)

    odesolver = lsrk

    vtkstep = [0, 0]
    cbvector =
        make_callbacks(vtkpath, vtkstep, nout, mpicomm, odesolver, dg, model, Q)

    eng0 = norm(Q)
    @info @sprintf """Starting
    norm(Q₀) = %.16e
    ArrayType = %s""" eng0 ArrayType

    solve!(Q, odesolver; timeend = timeend, callbacks = cbvector)

    Qe = init_ode_state(dg, timeend, init_on_cpu = true)

    error = euclidean_distance(Q, Qe) / norm(Qe)

    println("2D error = ", error)
    @test isapprox(error, FT(0.0); atol = 0.005)

    return nothing
end

function make_callbacks(
    vtkpath,
    vtkstep,
    nout,
    mpicomm,
    odesolver,
    dg,
    model,
    Q,
)
    if isdir(vtkpath)
        rm(vtkpath, recursive = true)
    end
    mkpath(vtkpath)
    mkpath(vtkpath * "/fast")

    function do_output(span, vtkstep, model, dg, Q)
        outprefix = @sprintf(
            "%s/%s/mpirank%04d_step%04d",
            vtkpath,
            span,
            MPI.Comm_rank(mpicomm),
            vtkstep
        )
        @info "doing VTK output" outprefix
        statenames = flattenednames(vars_state(model, Prognostic(), eltype(Q)))
        auxnames = flattenednames(vars_state(model, Auxiliary(), eltype(Q)))
        writevtk(outprefix, Q, dg, statenames, dg.state_auxiliary, auxnames)
    end

    do_output("fast", vtkstep[2], model, dg, Q)
    cbvtk = GenericCallbacks.EveryXSimulationSteps(nout) do (init = false)
        do_output("fast", vtkstep[2], model, dg, Q)
        vtkstep[2] += 1
        nothing
    end

    starttime = Ref(now())
    cbinfo = GenericCallbacks.EveryXWallTimeSeconds(60, mpicomm) do (s = false)
        if s
            starttime[] = now()
        else
            energy = norm(Q)
            @info @sprintf(
                """Update
                simtime = %8.2f / %8.2f
                runtime = %s
                norm(Q) = %.16e""",
                ODESolvers.gettime(odesolver),
                timeend,
                Dates.format(
                    convert(Dates.DateTime, Dates.now() - starttime[]),
                    Dates.dateformat"HH:MM:SS",
                ),
                energy
            )
        end
    end

    return (cbinfo)
end

#################
# RUN THE TESTS #
#################
FT = Float64
vtkpath =
    abspath(joinpath(ClimateMachine.Settings.output_dir, "vtk_bickley_jet"))

const timeend = FT(200) # s
const nout = FT(100)
const dt = FT(0.02) # s

const N = 3
const Nˣ = 8
const Nʸ = 8
const Lˣ = 4 / π  # m
const Lʸ = 4 / π  # m

xrange = range(FT(0); length = Nˣ + 1, stop = Lˣ)
yrange = range(FT(0); length = Nʸ + 1, stop = Lʸ)

run_bickley_jet(dt, nout)
