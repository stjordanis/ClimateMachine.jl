
using DiffEqOperators,  DifferentialEquations,  Plots
using Parameters
using Revise
#from jupyterthemes import jtplot
#jtplot.style(theme=’monokai’, context=’notebook’, ticks=True, grid=False)

FT = Float64

include("Water/water_functions.jl")
include("Water/hydraulic_conductivity_composable.jl")
include("Water/matric_potential_composable.jl")

hydraulic_model = vanGenuchten{FT}();
matric_model    = vanGenuchten{FT}();

K_sat  = 0.0443 / (3600*100)
porosity = 0.495 # Read in from data base
S_s = 10e-4  # [ m-1] ## Check this value !!

#IC/BC values  - to be specified/calculated via a BC/IC struct, and treated as an attribute of the model.
ν_0 = 0.24
ν_surface = porosity#+1e-5#-1e-3
S_l_0 = effective_saturation(porosity, ν_0)
ψ_0 = pressure_head(matric_model,S_l_0,porosity,S_s,ν_0)
println(ψ_0)
κ_0 = hydraulic_conductivity(hydraulic_model, K_sat, S_l_0, ψ_0,0.0)

L  = 1.0 # depth in meters
nknots = 99 #99
dz = -L/(nknots+1)#1000.0/(nknots+1)
knots = range(dz, step=dz, length=nknots)

z = collect(knots)
S = zeros(length(z))
#S[60]=0.0;#0.00001
ord_deriv = 2
ord_approx = 2

const Δ = CenteredDifference(ord_deriv, ord_approx, dz, nknots)
const ∇ = CenteredDifference(1, ord_approx, dz, nknots)
bc = DirichletBC(ν_surface, ν_0)


bc1 = RobinBC((1., 0., ν_surface), (0., 1., 0.), dz, 1)

function soil_water(ν,p,t; bc=bc, porosity=porosity,S=S, S_s=S_s, HM=hydraulic_model, MM=matric_model, K_sat=K_sat, z=z )
    # Get effective saturation
    ν_ext = bc*ν
    S_l = map(x -> effective_saturation(porosity,x),ν_ext)
    #@show S_l
    # This function calculates pressure head ψ of a soil
    ψ = similar(S_l)
    K = similar(S_l)
    for i in eachindex(ψ)
        ψ[i] = pressure_head(MM, S_l[i],porosity,S_s,ν_ext[i])
    end

    # Get hydraulic head
    #h = hydraulic_head(z,ψ)

    # Conductivity
    for i in eachindex(ψ)
        K[i] = hydraulic_conductivity(HM,K_sat,S_l[i],ψ[i],0.0)
    end
    #@show K
    dKdz = ∇*K

    K[2:end-1]*Δ*ψ + (∇*ψ) .* dKdz + dKdz + S

end

ν = ν_0*ones(99);
soil_water(ν,1,0);

tmax = 1*86400.0
prob = ODEProblem(soil_water, ν, (0.0, tmax))
alg = KenCarp4(autodiff=false)
sol = solve(prob, alg, saveat=10*60.0)


#sol_times = [sol(t) for t in 0.0:dt:tmax]

times = collect(0.:86400.0/4:tmax)
plt = plot()
for time in times
    plot!(plt, sol(time),z, xlabel="Theta_l", ylabel="Soil depths (m)", label=string("t(days) = ", time/86400.0), lw=3, legend=:bottomright)
end
for time in times
    plot!(plt, sol2(time),z, xlabel="Theta_l", ylabel="Soil depths (m)", label=string("t(days) = ", time/86400.0), lw=1, legend=:bottomright)
end
plot(plt, ylims=(-1.,0.))

#const bc = DirichletBC(porosity+1e-3, ν_0)
bc = DirichletBC(porosity+1e-4, ν_0)
@time prob = ODEProblem(soil_water, ν, (0.0, tmax))
alg = KenCarp4(autodiff=false)
sol2 = solve(prob, alg, saveat=10*60.0)

times = collect(0.:86400.0/5:tmax)
plt = plot()
for time in times
    plot!(plt, sol2(time),z, xlabel="Theta_l", ylabel="Soil depths (m)", label=string("t(days) = ", time/86400.0), lw=3, legend=:bottomright)
end

plot(plt, ylims=(-1.,0.))
