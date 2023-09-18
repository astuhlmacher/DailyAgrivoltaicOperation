# Select the project
cd(@__DIR__)

# Set environment
using Pkg
Pkg.activate(".")

# Load packages
using JuMP
using FileIO
using MAT
using Gurobi
using DataFrames
using CSV
using LinearAlgebra

# Import
include("src/apv_daily_problem.jl")

#################################################################
    # Input Data
#################################################################

## Time Step Information
π_e = repeat([0.03091, 0.02587, 0.02473, 0.02448, 0.02529, 0.02742, 0.02626, 0.02595, 0.02845, 0.03136, 0.03533, 0.03693, 0.04437, 0.06191, 0.06088, 0.05305, 0.04739, 0.04435, 0.04184, 0.03714, 0.03811, 0.03068, 0.02849, 0.02955], inner=(6,))  # electricity price over time [Units: $/kWh]
Δtpv = 10/60 # how frequently the PV position is updated [Units: h]

## PV Information
n_rows = 7 # number of rows [Units: -]
n_rowpairs = 20 # number of PV pairs within a row [Units: -]
row_distance = 6 #  distance between rows [Units: m]
pair_distance = 2.5 # distance between panel groups within row [Units: m]
PV_height = 4.5 # height of arrays above ground [Units: m]
PV_width = 1.135 # width of panel pair in E-W axis [Units: m]
PV_length = 4.2 # width of panel pair in N-S axis [Units: m]
η_array = 0.22 # PV efficiency [Units: -]
n_pvpairs = n_rows*n_rowpairs; # number of PV pairs [Units: -]
A_array = (n_pvpairs*PV_length*PV_width) # PV surface area [Units: m^2]

## Load Shading Analysis and Solar Position Information
filepath_sf = "data/AdjustTilt_90range_13density_fullday.mat";
sf_file = matopen(filepath_sf)
sun_azimuth = read(sf_file, "sun_azimuth")
sun_altitude = read(sf_file, "sun_altitude")
sf_coeffs = read(sf_file, "SF_coeffs") # Load linear shading factor deviation coefficients [a b] -> δ SF = ax+b, size 89x2
time_list = read(sf_file, "time_list")
daylight_start_index = 40 
daylight_end_index = 126
daylight_steps = length(daylight_start_index:daylight_end_index)
close(sf_file)

## Crop Information
crop_width = (PV_width + pair_distance)*n_rowpairs # [Units: meters]
crop_length = (PV_length+row_distance)*n_rows # [Units: meters]
PV_density = A_array/(crop_width*crop_length) # pv surface area / crop surface area [Units: -]
GHI2PAR = 0.44 # Ratio of PARtotal to GHI [Units: -]

crop_choice = 1 # Crop Choice 1: lettuce, 2: tomato
crop_LSP = [213.1, 596.6] # Light saturation point for various crops [Units: W/m^2]
daylight_hours_at_LSP =(time_list[daylight_end_index]-time_list[daylight_start_index]) # [Units: h]
PAR_req = 0.8*crop_LSP[crop_choice]*daylight_hours_at_LSP # Required PAR irradiance over the entire day [Units: Wh/m^2]

## Load Solar Radiance Information
filepath_irradiance = "data/AnnArbor_2021_simple_july14_EDT.csv"
irradiance_data = DataFrame(CSV.File(filepath_irradiance)) 
GHI = repeat(irradiance_data[:,:GHI], inner=(3,1)) # every 10 minutes, entire day [Units: W/m^2]
DHI = repeat(irradiance_data[:,:DHI], inner=(3,1)) # every 10 minutes, entire day [Units: W/m^2]
DNI = repeat(irradiance_data[:,:DNI], inner=(3,1)) # every 10 minutes, entire day [Units: W/m^2]

## Full Sun Tracking (ST) PV Power Output, for every 10 minutes in day 
PVenergy_ST,PVenergy_STtot = PVenergy_fullsuntrack(sun_altitude, sun_azimuth, GHI, DHI, DNI, Δtpv, A_array, η_array) # [Units: kWh]

## Full Sun Tracking (ST) daily accumulated PAR
PAR_ST = cropPAR_fullsuntrack(GHI, DHI, DNI, Δtpv, sf_coeffs, GHI2PAR) # [Units: Wh/m^2]

## Crop Daily Cumulative PAR in crop-only scenario
PAR_croponly = sum(GHI2PAR*GHI[t]*Δtpv for t=1:size(sun_altitude,1)) # [Units: Wh/m^2]


#################################################################
    # Solve Optimization
#################################################################

if !(@isdefined GRB_ENV)
    const GRB_ENV = Gurobi.Env()
end

model, fval, status = APV_optimization(sun_altitude, GHI, DHI, DNI, sf_coeffs, GHI2PAR, A_array, η_array, Δtpv, PAR_req, π_e, daylight_steps, daylight_start_index)

## Recover the Tilt Angle and Tilt Angle Deviation
x = value.(model[:x])
y = value.(model[:y])
test =  value.(model[:x]).^2 .+ value.(model[:y]).^2;
percent_inexact = sum(test.<0.999)/length(test)*100 # Evaluate exactness of recovered tilt angle values [%]
Nt = size(time_list,1)
δΣpv = zeros(Nt,1)
δΣpv[daylight_start_index:daylight_end_index] = sign.(asind.(y)).*maximum([abs.(acosd.(x)) abs.(asind.(y))], dims=2) # # Deviation of tilt angle from full sun tracking position
Σpv_ST = [min(max(i,0),90) for i in 90*ones(144,1)-sun_altitude] # Full sun-tracking tilt angle 
Σpv = [min(max(i,0),90) for i in 90*ones(144,1)-sun_altitude+δΣpv] # APV tilt angle within lower and upper tilt limits
Σpv_noclip = 90*ones(144,1)-sun_altitude+δΣpv # APV tilt angle with out lower and upper tilt limits clipping

## PV energy produced in optimal operation [Units: kWh]
PVenergy_opt = PVenergy_ST
PVenergy_opt[daylight_start_index:daylight_end_index] = PVenergy_ST[daylight_start_index:daylight_end_index]+ value.(model[:δP])*Δtpv*(1/1000) 
PVenergy_opttot = sum(PVenergy_opt)

## Daily cumulative PAR in optimal operation [Units: Wh/m^2/day]
PAR_opt = sum(value.(model[:PAR_crop])[t]*Δtpv for t=1:length(daylight_start_index:daylight_end_index))

PVprice_ST = dot(π_e, PVenergy_ST)
PVprice_opt = dot(π_e, PVenergy_opt) 

