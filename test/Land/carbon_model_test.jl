# First, we'll load our pre-requisites:
#  - load external packages:
using MPI
using OrderedCollections
using Plots
using StaticArrays

#  - load CLIMAParameters and set up to use it:

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

#  - load necessary ClimateMachine modules:
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils

#  - import necessary ClimateMachine modules: (`import`ing enables us to
#  provide implementations of these structs/methods)
import ClimateMachine.DGMethods:
    vars_state_auxiliary,
    vars_state_conservative,
    vars_state_gradient,
    vars_state_gradient_flux,
    source!,
    flux_second_order!,
    flux_first_order!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    update_auxiliary_state!,
    nodal_update_auxiliary_state!,
    init_state_auxiliary!,
    init_state_conservative!,
    boundary_state!

# ## Initialization

# Define the float type (`Float64` or `Float32`)
FT = Float64;
# Initialize ClimateMachine for CPU.
ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));

# Load some helper functions for plotting
include(joinpath(clima_dir, "docs", "plothelpers.jl"));

# # Define the set of Partial Differential Equations (PDEs)

# ## Define the model

# Model parameters can be stored in the particular [`BalanceLaw`](@ref
# ClimateMachine.DGMethods.BalanceLaw), in this case, a `CarbonModel`:

"""
    CarbonModel

Simple carbon model
"""
Base.@kwdef struct CarbonModel{FT} <: BalanceLaw
    "Parameters"
    param_set::AbstractParameterSet = param_set
    "Initial B (biomass) (g carbon/m^2)"
    B_init::FT = 5000
    "Initial S (soil) (g carbon/m^2)"
    S_init::FT = 20000
    "Net primary production (g carbon/m^2/yr)"
    NPP::FT = 295.15
    "k_1 (1/yr)"
    k_1::FT = 0.1*1
    "k_2 (1/yr)"
    k_2::FT = 0.02*1
end

# Create an instance of the `CarbonModel`:
m = CarbonModel{FT}();

# This model dictates the flow control, using [Dynamic Multiple
# Dispatch](https://en.wikipedia.org/wiki/Multiple_dispatch), for which
# kernels are executed.

# ## Define the variables

# All of the methods defined in this section were `import`ed in # [Loading
# code](@ref) to let us provide implementations for our `CarbonModel` as they
# will be used by the solver.

# Specify auxiliary variables for `CarbonModel`
vars_state_auxiliary(::CarbonModel, FT) = @vars(Sk_2::FT);

# Specify state variables, the variables solved for in the PDEs, for
# `CarbonModel`
vars_state_conservative(m::CarbonModel, FT) = @vars(B::FT, S::FT);

# Specify state variables whose gradients are needed for `CarbonModel`
vars_state_gradient(::CarbonModel, FT) = @vars();

# Specify gradient variables for `CarbonModel`
vars_state_gradient_flux(::CarbonModel, FT) = @vars();

# ## Define the compute kernels

# Specify the initial values in `aux::Vars`, which are available in
# `init_state_conservative!`. Note that
# - this method is only called at `t=0`
# - `aux.z` and `aux.T` are available here because we've specified `z` and `T`
# in `vars_state_auxiliary`
function init_state_auxiliary!(m::CarbonModel, aux::Vars, geom::LocalGeometry)
    aux.Sk_2 = m.S_init*m.k_2
end;

# Specify the initial values in `state::Vars`. Note that
# - this method is only called at `t=0`
# - `state.ρcT` is available here because we've specified `ρcT` in
# `vars_state_conservative`
function init_state_conservative!(
    m::CarbonModel,
    state::Vars,
    aux::Vars,
    coords,
    t::Real,
)
    state.B = m.B_init
    state.S = m.S_init
end;

# The remaining methods, defined in this section, are called at every
# time-step in the solver by the [`BalanceLaw`](@ref
# ClimateMachine.DGMethods.BalanceLaw) framework.

# Overload `update_auxiliary_state!` to call `carbon_eq_nodal_update_aux!`, or
# any other auxiliary methods
function update_auxiliary_state!(
    dg::DGModel,
    m::CarbonModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
    nodal_update_auxiliary_state!(carbon_eq_nodal_update_aux!, dg, m, Q, t, elems)
    return true # TODO: remove return true
end;

# Compute/update all auxiliary variables at each node. Note that
# - `aux.T` is available here because we've specified `T` in
# `vars_state_auxiliary`
function carbon_eq_nodal_update_aux!(
    m::CarbonModel,
    state::Vars,
    aux::Vars,
    t::Real,
)
    aux.Sk_2 = state.S * m.k_2
end;

# Since we have second-order fluxes, we must tell `ClimateMachine` to compute
# the gradient of `ρcT`. Here, we specify how `ρcT` is computed. Note that
#  - `transform.ρcT` is available here because we've specified `ρcT` in
#  `vars_state_gradient`
function compute_gradient_argument!(
    m::CarbonModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
end;

# Specify where in `diffusive::Vars` to store the computed gradient from
# `compute_gradient_argument!`. Note that:
#  - `diffusive.α∇ρcT` is available here because we've specified `α∇ρcT` in
#  `vars_state_gradient_flux`
#  - `∇transform.ρcT` is available here because we've specified `ρcT`  in
#  `vars_state_gradient`
function compute_gradient_flux!(
    m::CarbonModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
end;

# We have no sources, nor non-diffusive fluxes.
function source!(
    m::CarbonModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    Bk_1 = state.B*m.k_1
    source.B = m.NPP - Bk_1
    source.S = Bk_1 - aux.Sk_2
end;
function flux_first_order!(
    m::CarbonModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
) end;

# Compute diffusive flux (``F(α, ρcT, t) = -α ∇ρcT`` in the original PDE).
# Note that:
# - `diffusive.α∇ρcT` is available here because we've specified `α∇ρcT` in
# `vars_state_gradient_flux`
function flux_second_order!(
    m::CarbonModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
end;

# ### Boundary conditions

# Second-order terms in our equations, ``∇⋅(F)`` where ``F = -α∇ρcT``, are
# internally reformulated to first-order unknowns.
# Boundary conditions must be specified for all unknowns, both first-order and
# second-order unknowns which have been reformulated.

# The boundary conditions for `ρcT` (first order unknown)
function boundary_state!(
    nf,
    m::CarbonModel,
    state⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
end;

# The boundary conditions for `ρcT` are specified here for second-order
# unknowns
function boundary_state!(
    nf,
    m::CarbonModel,
    state⁺::Vars,
    diff⁺::Vars,
    aux⁺::Vars,
    n⁻,
    state⁻::Vars,
    diff⁻::Vars,
    aux⁻::Vars,
    bctype,
    t,
    _...,
)
end;

# # Spatial discretization

# Prescribe polynomial order of basis functions in finite elements
N_poly = 5;

# Specify the number of vertical elements
nelem_vert = 3;

# Specify the domain height
zmax = FT(1);

# Establish a `ClimateMachine` single stack configuration
driver_config = ClimateMachine.SingleStackConfiguration(
    "CarbonEquation",
    N_poly,
    nelem_vert,
    zmax,
    param_set,
    m,
    numerical_flux_first_order = CentralNumericalFluxFirstOrder(),
);

# # Time discretization

# Specify simulation time (SI units)
t0 = FT(0)
timeend = FT(100)
dt = FT(1)

# # Configure a `ClimateMachine` solver.

# This initializes the state vector and allocates memory for the solution in
# space (`dg` has the model `m`, which describes the PDEs as well as the
# function used for initialization). This additionally initializes the ODE
# solver, by default an explicit Low-Storage
# [Runge-Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
# method.

solver_config =
    ClimateMachine.SolverConfiguration(t0, timeend, driver_config, ode_dt = dt);
grid = solver_config.dg.grid;
Q = solver_config.Q;
aux = solver_config.dg.state_auxiliary;

# ## Inspect the initial conditions

# Let's export a plot of the initial state
output_dir = @__DIR__;

mkpath(output_dir);

z_scale = 100 # convert from meters to cm
z_key = "z"
z_label = "z [cm]"
z = get_z(grid, z_scale)
state_vars = SingleStackUtils.get_vars_from_nodal_stack(
    grid,
    Q,
    vars_state_conservative(m, FT),
)
aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
    grid,
    aux,
    vars_state_auxiliary(m, FT),
)
all_vars = OrderedDict(state_vars..., aux_vars...);
export_plot_snapshot(
    z,
    all_vars,
    ("B", "S",),
    joinpath(output_dir, "initial_condition.png"),
    z_label,
);
# ![](initial_condition.png)

# It matches what we have in `init_state_conservative!(m::CarbonModel, ...)`, so
# let's continue.

# Solver hooks / callbacks

# Define the number of outputs from `t0` to `timeend`
const n_outputs = 20;

# This equates to exports every ceil(Int, timeend/n_outputs) time-step:
const every_x_simulation_time = ceil(Int, timeend / n_outputs);

# Create a nested dictionary to store the solution:
all_data = OrderedDict([k => Dict() for k in 1:n_outputs+1]...)
t_vec = FT[0 for i in 1:n_outputs+1]

# The `ClimateMachine`'s time-steppers provide hooks, or callbacks, which
# allow users to inject code to be executed at specified intervals. In this
# callback, the state and aux variables are collected, combined into a single
# `OrderedDict` and written to a NetCDF file (for each output step `step`).
step = [1];
callback = GenericCallbacks.EveryXSimulationTime(
    every_x_simulation_time,
    solver_config.solver,
) do (init = false)
    state_vars = SingleStackUtils.get_vars_from_nodal_stack(
        grid,
        Q,
        vars_state_conservative(m, FT),
    )
    aux_vars = SingleStackUtils.get_vars_from_nodal_stack(
        grid,
        aux,
        vars_state_auxiliary(m, FT);
        exclude = [z_key],
    )
    all_vars = OrderedDict(state_vars..., aux_vars...)
    all_data[step[1]] = all_vars
    t_vec[step[1]] = ODESolvers.gettime(solver_config.solver)

    step[1] += 1
    nothing
end;

# # Solve

# This is the main `ClimateMachine` solver invocation. While users do not have
# access to the time-stepping loop, code may be injected via `user_callbacks`,
# which is a `Tuple` of [`GenericCallbacks`](@ref).
ClimateMachine.invoke!(solver_config; user_callbacks = (callback,));

# # Post-processing

S_vs_t = [all_data[i]["S"][1] for i in keys(all_data)]
B_vs_t = [all_data[i]["B"][1] for i in keys(all_data)]
plot(t_vec, B_vs_t, label = "B")
plot!(t_vec, S_vs_t, label = "S")
savefig(joinpath(output_dir, "sol_vs_time.png"))
# File is now saved in `output_dir`
