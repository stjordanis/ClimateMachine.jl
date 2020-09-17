const clima_dir = dirname(dirname(pathof(ClimateMachine)));
output_dir = joinpath(clima_dir, "output")
gcm_output_dir = joinpath(output_dir, "gcm")
fig_dir = joinpath(gcm_output_dir, "figs")
mkpath(fig_dir)

all_files = readdir(output_dir)
all_files = map(x->joinpath(output_dir, x), all_files)
nc_files = filter(f->endswith(f, ".nc"), all_files)
if isempty(nc_files)
    @warn "Using the same data"
else
    @info "new data being used"
    data_file = nc_files[1]
    mv(data_file, joinpath(gcm_output_dir,"gcm_output.nc"); force=true)
    data_file = joinpath(gcm_output_dir,"gcm_output.nc")
    if isfile(data_file)
        for f in nc_files
            rm(f; force=true)
        end
    else
        error("Oops")
    end
end

using NCDatasets
ds = Dataset(data_file, "r")

function get_column(ds, var_name, i_time, i_lat, i_long)
    lat = ds["lat"][:][i_lat]
    long = ds["long"][:][i_long]
    z = ds["level"][:] / 1e3
    t = ds["time"][:][i_time]
    var = ds[var_name][:]
    data = var[i_long,i_lat,:,i_time] # long, lat, lev, time
    return z, lat, long, data, t
end
using Plots

function get_plots(ds, v_all, i_time, fig_dir, plot_single=true)
    plot_array = Any[];
    for v in v_all
        title = v
        one_plot = plot()
        for i_lat in 1:30:ds.dim["lat"]
            for i_long in 1:60:ds.dim["long"]
                z, lat, lon, data, t = get_column(ds, v, i_time, i_lat, i_long)
                one_plot = plot!(data, z, title = v, ylabel="height (km)", xlabel=v, label="lat/lon=$(lat),$(lon)", legend=:topleft);
            end
        end
        plot_single && savefig(string(joinpath(fig_dir,"$v vs z.pdf")))
        push!(plot_array,one_plot); # make a plot and add it to the plot_array
    end
    return plot_array
end

i_time = 6
v_all = ["rho", "temp", "pres"]
plot_array = get_plots(ds, v_all, i_time, fig_dir, true)

fig=plot(plot_array... , layout=(1,length(v_all)), size=(1000, 400) )
savefig(fig, string(joinpath(fig_dir,"ref_summary.pdf")))

v_all = ["thd", "et", "ei", "ht", "hi"]
plot_array = get_plots(ds, v_all, i_time, fig_dir, false)
fig=plot(plot_array... , layout=(1,length(v_all)), size=(1000, 400) )
savefig(fig, string(joinpath(fig_dir,"energy_summary.pdf")))

v_all = ["u","v", "w"]
plot_array = get_plots(ds, v_all, i_time, fig_dir, true)
