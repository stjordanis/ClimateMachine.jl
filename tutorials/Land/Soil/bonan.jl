# ---------------------------------------------------------------------
# Use the predictor-corrector method to solve the Richards equation for
# infiltration with surface soil moisture as the boundary condition.
# ---------------------------------------------------------------------
# Set The PATH Here!
#push!(LOAD_PATH, "/Users/ddutta/Dropbox/JULIA_Research/Soil_Moisture_Codes/PREDICTOR_CORRECTOR/RichardsEquation/")

using Plots
using Parameters
using soil_moisture_V4
global soil = soil_moisture_V4.soil_struct(nsoi=150) # Define the number of soil Layers

# Generate the Soil Compute Grid
soil = soil_moisture_V4.compute_grid_settings(soil)

# --- Soil parameters Settings
soil.ssflag    = 0 # Flag for Sink term
soil.functions = "van_Genuchten";  # Use van Genuchten relationships
#soil.functions = "Campbell";       # Use Campbell relationships

if soil.functions == "Campbell"

   # example from Hornberger & Wiberg [2005, Fig. 8.3]
   ityp = 0;              # Soil texture flag
   theta_sat = 0.25;      # Volumetric water content at saturation
   psi_sat = -25.0;       # Matric potential at saturation [cm]
   bc = 0.2;              # Exponent
   Ksat = 3.4e-03;        # Hydraulic conductivity at saturation [cm/s]
   params = [theta_sat psi_sat bc Ksat ityp]

elseif soil.functions == "van_Genuchten"

   # Haverkamp et al. (1977): sand
   ityp = 1;              # Soil texture flag
   theta_res = 0.075;     # Residual water content
   theta_sat = 0.287;     # Volumetric water content at saturation
   vg_alpha = 0.027;      # Inverse of the air entry potential [/cm]
   vg_n = 3.96;           # Pore-size distribution index
   vg_m = 1;              # Exponent
   Ksat = 34 / 3600;      # Hydraulic conductivity at saturation [cm/s]

#  # Haverkamp et al. (1977): Yolo light clay
#  ityp = 2;              # Soil texture flag
#  theta_res = 0.124;     # Residual water content
#  theta_sat = 0.495;     # Volumetric water content at saturation
#  vg_alpha = 0.026;      # Inverse of the air entry potential [/cm]
#  vg_n = 1.43;           # Pore-size distribution index
#  vg_m = 1 - 1 / vg_n;   # Exponent
#  Ksat = 0.0443 / 3600;  # Hydraulic conductivity at saturation [cm/s]

   params = [theta_res theta_sat vg_alpha vg_n vg_m Ksat ityp]
end

# --- Initial soil moisture & matric potential
for i = 1:soil.nsoi
   if (ityp == 0)
      soil.theta[i] = 0.10
   elseif ityp == 1
      soil.theta[i] = 0.10
   elseif ityp == 2
      soil.theta[i] = 0.24
   end
   soil.psi[i] = matric_potential(soil.functions, params, soil.theta[i])
end

global thetaini = copy(soil.theta)

# --- Surface boundary condition: saturation [minus some small delta]

soil.theta0 = theta_sat - 1.0e-03
if (ityp == 1)
   soil.theta0 = 0.267
end
soil.psi0 = matric_potential(soil.functions, params, soil.theta0)

# --- Time step [seconds]

dt = 10
if (ityp == 1)
   dt = 5
end

# --- Length of simulation [number of time steps]

# Hornberger & Wiberg: 15; 30; | 60 minutes
if (ityp == 0)
#  ntim = 15 * 60 / dt
#  ntim = 30 * 60 / dt
   ntim = 60 * 60 / dt
end

# Haverkamp et al. (1977) - sand: duration is in hours
if (ityp == 1)
#  ntim = 0.05 * 3600 / dt
#  ntim = 0.1 * 3600 / dt
#  ntim = 0.2 * 3600 / dt
#  ntim = 0.3 * 3600 / dt
# ntim = 0.4 * 3600 / dt
  ntim = 0.8 * 3600 / dt
end

# Haverkamp et al. (1977) - Yolo light clay: duration is in seconds
if (ityp == 2)
#  ntim = 1.0e4 / dt
#  ntim = 1.0e5 / dt
#  ntim = 5.0e5 / dt
   ntim = 1.0e6 / dt
end


# --- Initialize accumulators for water balance check
global sum_in    = 0
global sum_out   = 0
global sum_store = 0

global ET = 0.006

# --- Time stepping loop: NTIM iterations with a time step of DT seconds
# Initialize cumulative infiltration variables
xout = zeros(Int64.(ntim))
yout = zeros(Int64.(ntim))
for itim = 1:Int64.(ntim)

   # Hour of day

   hour = itim * (dt/86400 * 24)
   @printf("hour = %8.3f \n",hour)
    
#  # Add for fun variable Boundary Condition SURFACE DRYING
#     if itim > 100
#         soil.theta0 = theta_sat - 1.0e-03
#         if (ityp == 1)
#            soil.theta0 = 0.267 - (itim-200)*0.1/ntim
#         end
#         soil.psi0 = matric_potential(soil.functions, params, soil.theta0)
#     end


   # Calculate soil moisture
   soil = predictor_corrector(soil, params, ET, dt)


   # Sum fluxes for relative mass balance error()

   sum_in = sum_in + abs(soil.Q0) * dt
   sum_out = sum_out + abs(soil.QN) * dt
   sum_store = sum_store + soil.dtheta
    if itim%10 == 0
        IJulia.clear_output(true)
        Plots.display(plot(soil.theta,soil.z, color="red",line=2,xlabel = "theta [-]",ylabel = "z [cm]",title="Time = "*string(round(hour*60,digits=2))*" mins"))
        sleep(0.2)
    end
        
    #print(soil.sink)

   # cumulative infiltration
   xout[itim] = hour
   yout[itim] = sum_in

end