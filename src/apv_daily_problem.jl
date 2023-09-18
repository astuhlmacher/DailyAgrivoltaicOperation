
function PVenergy_fullsuntrack(sun_altitude, sun_azimuth, GHI, DHI, DNI, Δtpv, A_array, η_array)
    ## Assuming no shading between panels, which may be very reasonable since agrivoltaic panels use less dense configurations
    ## Full Sun Tracking 
    pv_azimuth_track = sun_azimuth # PV azimuth position [degrees]
    pv_tilt_track = 90 .- sun_altitude # PV tilt position [degrees] 
    kilo = 1000 # scaling factor

    PVenergy_ST = zeros(size(sun_altitude,1),1) # init PV energy value
    for t=1:size(sun_altitude,1)

        IBt = DNI[t]  # direct beam irradiance, zero incidence angle (changes every 30 minutes, save data every 10 minutes)
        IDt = DHI[t]*((1+cosd(pv_tilt_track[t]))/2)# diffuse beam irradiance (pv updates every 10 minutes)

        PVenergy_ST[t] = A_array*η_array*(IBt + IDt)*Δtpv/kilo # [kWh]
    end

    PVenergy_STtot = sum(PVenergy_ST)

    return PVenergy_ST, PVenergy_STtot
end

function cropPAR_fullsuntrack(GHI, DHI, DNI, Δtpv, sf_coeffs, GHI2PAR)
    ## Full Sun Tracking 
    #pv_azimuth_track = sun_azimuth # PV azimuth position [degrees]
    #pv_tilt_track = 90 .- sun_altitude # PV tilt position [degrees] 
    # δϕ = 0 -> cos(δϕ) =1 -> x=1

    cropPAR_ST = 0 # init cumulative PAR value
    for t=1:size(sf_coeffs,1) 
        PAR_t = GHI2PAR*GHI[t]-sf_coeffs[t,1]*GHI2PAR*GHI[t]-sf_coeffs[t,2]*GHI2PAR*GHI[t]
        
        cropPAR_ST += PAR_t*Δtpv # [Wh/m^2]
    end

    return cropPAR_ST
end

function cropSF(sf_coeffs, Nt, δΣpv)
    ## optimal Tracking 

    cropSF_opt = zeros(Nt,1) # init cumulative PAR value
    for t=1:Nt
        cropSF_opt[t] = sf_coeffs[t,1]*cosd(δΣpv[t])+sf_coeffs[t,2]
    end
    return cropSF_opt
end


function APV_optimization(sun_altitude, GHI, DHI, DNI, sf_coeffs, GHI2PAR, A_array, η_array, Δtpv, PAR_req, π_e, Nt, daylight_start_index)
    # Init Model with Solver
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "LogtoConsole", 0)


    ## Init Decision Variables
    @variable(model, x[t=1:Nt], lower_bound=0, upper_bound=1, base_name="x_tilt")
    @variable(model, y[t=1:Nt], lower_bound=-1, upper_bound=1, base_name="y_tilt") 
    @variable(model, δIB[t=1:Nt], base_name="dev IB")
    @variable(model, δID[t=1:Nt], base_name="dev ID")
    @variable(model, δP[t=1:Nt], base_name="dev Power")
    @variable(model, PAR_crop[t=1:Nt],  base_name="PARcrop")

    ## Define a, b, c, and d parameter matrices
    a = zeros(2,Nt)
    b = zeros(3,Nt)
    c = zeros(2,Nt)
    d = zeros(2,Nt)
    f = zeros(2,Nt)
    for t=1:Nt
        a[1,t] = DNI[daylight_start_index-1+t] 
        a[2,t] = -DNI[daylight_start_index-1+t] 
        b[1,t] = (1/2)*DHI[daylight_start_index-1+t]*cosd(90-sun_altitude[daylight_start_index-1+t]) 
        b[2,t] = -(1/2)*DHI[daylight_start_index-1+t]sind(90-sun_altitude[daylight_start_index-1+t]) 
        b[3,t] = -(1/2)*DHI[daylight_start_index-1+t]*cosd(90-sun_altitude[daylight_start_index-1+t]) 
        c[1,t] = cosd(90-sun_altitude[daylight_start_index-1+t]) 
        c[2,t] = -sind(90-sun_altitude[daylight_start_index-1+t])
        d[1,t] = -sf_coeffs[daylight_start_index-1+t,1]*GHI2PAR*GHI[daylight_start_index-1+t] 
        d[2,t] = GHI2PAR*GHI[daylight_start_index-1+t] - sf_coeffs[daylight_start_index-1+t,2]*GHI2PAR*GHI[daylight_start_index-1+t] 
        f[1,t] = sind(90-sun_altitude[daylight_start_index-1+t]) 
        f[2,t] = cosd(90-sun_altitude[daylight_start_index-1+t])
    end

    ## Power Deviations
    @constraint(model, PowerDev[t=1:Nt], δP[t] == A_array*η_array*(δIB[t] + δID[t]))

    ## Irradiance Deviations
    @constraint(model, IbeamDev[t=1:Nt], δIB[t] == a[1,t]*x[t] + a[2,t] )
    @constraint(model, IdiffDev[t=1:Nt], δID[t] == b[1,t]*x[t] + b[2,t]*y[t] + b[3,t] )

    ## x,y limits
    @constraint(model, XYlim_LB2[t=1:Nt], f[1,t]*x[t]+f[2,t]*y[t] >= 0)
    @constraint(model, XYlim_UB2[t=1:Nt], f[1,t]*x[t]+f[2,t]*y[t] <= 1)
    @constraint(model, XY_relationship[t=1:Nt], x[t]^2 + y[t]^2 <=1)

    ## PAR definition and requirements
    @constraint(model, PARcrop_def[t=1:Nt], PAR_crop[t] == d[1,t]*x[t] + d[2,t])
    @constraint(model, PAR_req, sum(PAR_crop[t]*Δtpv for t=1:Nt) >= PAR_req)

    ## Objective Function
    @objective(model, Max, sum(π_e[daylight_start_index-1+t]*δP[t]*Δtpv  for t=1:Nt) )

    optimize!(model) # Solves optimization problem
    @show status = termination_status(model) # Determines the solver's status
    
    fval=NaN
    try
        fval=objective_value(model)
    catch
        println(string("Error: ",status))
        fval=Inf
    end

    return model, fval, status 
end