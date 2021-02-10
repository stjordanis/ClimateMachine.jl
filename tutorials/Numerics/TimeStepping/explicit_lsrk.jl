# # Single-rate explicit timestepping

# ClimateMachine supports multiple timesteppers, of different nature (explicit, semi-implicit, single-stage, multi-stage, single-step, multi-step, single-rate, multi-rate, etc). In this tutorial, we shall explore the use of explicit Runge-Kutta methods. For our model problem, we shall reuse the [rising thermal bubble tutorial](@ref Rising-termal-bubble)
# (see for details on the model and parameters)


include("tutorials/Numerics/TimeStepping/tutorial_risingbubble_config.jl")
FT = Float64

# After discretizing in space, the semi-discretization of the governing equations have the form:

# $$
# \begin{aligned}
#     \frac{\mathrm{d} \boldsymbol{q}}{ \mathrm{d} t} &= M^{-1}\left(M S +
#     D^{T} M (F^{adv} + F^{visc}) + \sum_{f=1}^{N_f} L^T M_f(\widehat{F}^{adv} + \widehat{F}^{visc})
#     \right) \equiv \mathcal{T}(\boldsymbol{q}).
# \end{aligned}
# $$
#
# # ## Low-storage Runge-Kutta methods
#
# [Intro something along these lines]
#
# A single step of an ``s``-stage Runge-Kutta (RK) method for
# solving the resulting ODE problem in [eq:foo] and can be
# expressed as the following:
#
# $$
# \begin{aligned}
# 	\boldsymbol{q}^{n+1} = \boldsymbol{q}^n + \Delta t \sum_{i=1}^{s} b_i \mathcal{T}(\boldsymbol{Q}^i),
# \end{aligned}
# $$
#
# where $\boldsymbol{\mathcal{T}}(\boldsymbol{Q}^i)$ is the evaluation of the right-hand side tendency at the stage value $\boldsymbol{Q}^i$, defined at each stage of the RK method:
#
# $$
# \begin{aligned}
# 	\boldsymbol{Q}^i = \boldsymbol{q}^{n} +
#     \Delta t \sum_{j=1}^{s} a_{i,j}
#     \mathcal{T}(\boldsymbol{Q^j}).
# \end{aligned}
# $$
#
# The first stage is initialized using the field at the previous time-step: $\boldsymbol{Q}^{1} = \boldsymbol{q}^n$.
#
# In the above expressions, we define $\boldsymbol{A} = \lbrace a_{i,j} \rbrace \in \mathbb{R}^{s\times s}$, $\boldsymbol{b} = \lbrace b_i \rbrace \in \mathbb{R}^s$, and $\boldsymbol{c} = \lbrace c_i \rbrace \in \mathbb{R}^s$ as the characteristic coefficients of a given RK method. This means we can associate any RK method with its so-called *Butcher tableau*:
#
# $$
# \begin{aligned}
#     \begin{array}{c|c}
#         \boldsymbol{c} &\boldsymbol{A}\\
#         \hline
#         & \boldsymbol{b}^T
#         \end{array} =
#         \begin{array}{c|c c c c}
#         c_1 & a_{1,1} & a_{1,2} & \cdots & a_{1,s}\\
#         c_2 & a_{2,1} & a_{2,2} & \cdots & a_{2,s}\\
#         \vdots & \vdots & \vdots & \ddots & \vdots\\
#         c_s & a_{s,1} & a_{s,2} & \cdots & a_{s,s}\\
#         \hline
#         & b_1 & b_2 & \cdots & b_s
#     \end{array}.
# \end{aligned}
# $$
#
# The vector $\boldsymbol{c}$ is often called the *consistency vector*, and is typically subject to the row-sum condition:
#
# $$
# c_i = \sum_{j=1}^{s} a_{i,j}, \quad \forall i = 1, \cdots, s.
# $$
#
# This simplifies the order conditions for higher-order RK methods. For more information on general RK methods, we refer the interested reader to \cite[Ch. 5.2]{atkinson2011numerical}.
#
# ## [Low-storage Runge-Kutta methods](@id lsrk)
#
# Here, we use a 4-th order 5-stage explict LSRK method:

ode_solver = ClimateMachine.ExplicitSolverType(
    solver_method = LSRK54CarpenterKennedy,
)
CFL = FT(0.6)
timeend = FT(100)
run_simulation(ode_solver, CFL, timeend)

# Let's try to take a larger time-step:

CFL = FT(1.7)
run_simulation(ode_solver, CFL, timeend)

# Oh-no it breaks!
# Now we try using a 14-stage method. Due to its larger
# stability region, we can take a larger time-step size
# compared to the previous 5-stage method:

ode_solver = ClimateMachine.ExplicitSolverType(
    solver_method = LSRK144NiegemannDiehlBusch,
)
run_simulation(ode_solver, CFL, timeend)

# # Strong Stability Preserving Runge--Kutta methods

# Just as with the LSRK methods, the SSPRK methods are self-starting, with $\boldsymbol{Q}^{1} = \boldsymbol{q}^n$, and stage-values are of the form

# $$
# \begin{aligned}
#     \boldsymbol{Q}^{i+1} = a_{i,1} \boldsymbol{q}^n
#     + a_{i,2} \boldsymbol{Q}^{i}
#     + \Delta t b_i\mathcal{T}(\boldsymbol{Q}^{i})
# \end{aligned}
# $$
#
# with the value at the next step being the $(N+1)$-th stage value $\boldsymbol{q}^{n+1} = \boldsymbol{Q}^{(N+1)}$. This allows the updates to be performed with only three copies of the state vector (storing $\boldsymbol{q}^n$, $\boldsymbol{Q}^{i}$ and $\mathcal{T}(\boldsymbol{Q}^{i})$).

# [Reference literature for the theoretical construction of the SSPRK methods]

# ## [Strong-stability-preserving Runge-Kutta methods](@id ssprk)

ode_solver = ClimateMachine.ExplicitSolverType(
    solver_method = SSPRK33ShuOsher,
)
CFL = FT(1.7)
run_simulation(ode_solver, CFL, timeend)
