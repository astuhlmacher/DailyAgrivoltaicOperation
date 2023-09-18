%%  Shading Calculations
% ADJUSTING TILT
clear all;

%% Operation Choice
tilt_limits = [0,90]; % degrees
azimuth_limits = [-180, 180]; % degrees
discrete_step = 1; % degrees
discrete_range = 90; % degrees

%% Sun Position
% Load solar position for specific day and year 
filepath = "july-14-2021data.mat";
load(filepath)

% Get solar position wrt definitions given in Cascone et al.
sun_azimuth = 180-sun_azimuth; % solar azimuth angle, angle of sun east of due south [degrees]...positive in morning, negative in afternoon 
sun_altitude = 90-sun_zenith; % solar altitude angle, angle of sun above horizon, complement of zenith angle [degrees]...if negative then you cannot see the sun from the crop area (ignore negative values)

%% PV Panels
% There are X rows of PV panels and Y pairs of PV panels per row.
% We define the location of the PV pairs from the center of each pair, where the height equals stilt or main axis height
% We define the panel decisions in terms of tilt and azimuth.
n_rows = 7; % number of rows
n_pvpairs = 20; % number of PV pairs within a row
row_distance = 6; %  distance between rows [m]
pair_distance = 2.5; % distance between panel groups within row [m]
PV_width = 1.135; % width in E-W axis[m]
PV_length = 4.2; % width in N-S axis [m]
PV_height = 4.5; % height of panels[m]
n_pv = n_rows*n_pvpairs; % number of PV pairs

PV_pairs = zeros(n_pv,3); % (X,Y)=(E,N)...coordinates changed later to SEZ global coordinate system
for row = 1:n_rows
    for pair = 1:n_pvpairs
        PV_pairs(n_pvpairs*(row-1)+pair,:) = [(pair_distance+PV_width)*((2*pair-1)/2), (row_distance+PV_length)*((2*row-1)/2), PV_height];
    end
end

%% Crop Land Info
% Crop land: SW vertex is (0,0) on local coordinate system
% Crop land has width (x) of XX m and length (y) of YY m
crop_width = (PV_width + pair_distance)*n_pvpairs; % meters
crop_length = (PV_length+row_distance)*n_rows; % meters
pv_density = (n_pv*PV_length*PV_width)/(crop_width*crop_length); % pv surface area / crop surface area

crop_azimuth = 0; % degrees

%% Time-varying Shading Factor
SF = {}; % Init time-varying shading factor, {t, 1:[tilt, azimuth], 2:[SFbase], 3:[tilt deviations],  4:[SF actual], 5:[SF_dev]}
SF_coeffs = zeros(size(sun_altitude,1),2);
for t=1:size(time_list,1) % Iterate over time horizon
    t  
    if sun_altitude(t) >= 0  % only calculate shading factor when sun is above horizon
        PV_tilt_des = max(min(90 - sun_altitude(t),tilt_limits(2)), tilt_limits(1));  % degrees, with clipping
        PV_azimuth_des = max(min(sun_azimuth(t), azimuth_limits(2)), azimuth_limits(1));  % degrees, with clipping
        SF{t,1} = [PV_tilt_des, PV_azimuth_des];
  
        SF_base = PV_shading_factor(sun_azimuth(t), sun_altitude(t), PV_azimuth_des, PV_tilt_des, n_pv, PV_pairs, PV_length, PV_width, crop_width, crop_length, crop_azimuth);      
        SF{t,2} = SF_base;
     
        % Range of Deviations Considered
        % Tilt Angle Deviations
        SF{t,3} = max(tilt_limits(1) - PV_tilt_des, -discrete_range):discrete_step:min(tilt_limits(2) - PV_tilt_des, discrete_range); 

        % Calculate shading factor changes given panel operation deviations
        SF{t,4} = []; %init  SF actual
        SF{t,5} = []; %init SF dev
        for ii = 1:size(SF{t,3},2) % go through all tilt angle deviations
            % only calculate shading factor when sun is above horizon
            SF{t,4} = [SF{t,4}, PV_shading_factor(sun_azimuth(t), sun_altitude(t), PV_azimuth_des, PV_tilt_des+SF{t,3}(ii), n_pv, PV_pairs, PV_length, PV_width, crop_width, crop_length, crop_azimuth)];      
            SF{t,5} = [SF{t,5}, SF{t,4}(ii)-SF_base];   
        end
        SF_coeffs(t,:) = polyfit(cosd(SF{t,3}),SF{t,4},1);
    else
        SF{t,1} = nan;
        SF{t,2} = nan;
        SF{t,3} = nan;
        SF{t,4} = nan;
        SF{t,5} = nan;
        SF_coeffs(t,1:2) = [0,0];
    end
end

%% Plot shading factor over optimization horizon
figure(1)
clf 
hold on;

for t=1:size(SF,1)
    plot(cosd(SF{t,3}), SF{t,5}, 'DisplayName', string(duration(hours(time_list(t)),'format','hh:mm')))
end
xlabel('Cosine Tilt Deviations from Sun Tracking (degrees)')
ylabel('Shading Factor Deviation')

figure(2)
clf 
hold on;

for t=1:size(SF,1)
    plot(cosd(SF{t,3}), SF{t,4}, 'DisplayName', string(duration(hours(time_list(t)),'format','hh:mm')))
end
xlabel('Cosine Tilt Deviations from Sun Tracking (degrees)')
ylabel('Shading Factor')

figure(3)
clf 
hold on; 

for ii=40:126; %1:size(SF,1)
    x_range = cosd(SF{ii,3}); %cosd(SF{ii,3}(1)):0.05:cosd(SF{ii,3}(end));
    plot(x_range, SF_coeffs(ii,1)*x_range + SF_coeffs(ii,2), 'DisplayName', string(duration(hours(time_list(ii)),'format','hh:mm')))
end
xlabel("cos tilt")
ylabel("linear best fit")
    
save("AdjustTilt_90range_13density_fullday.mat","SF","time_list","sun_altitude","sun_azimuth", "SF_coeffs", "time_list")

function SF = PV_shading_factor(sun_azimuth, sun_altitude, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, crop_width, crop_length, crop_azimuth)
    %% Outputs
    % SF: shading factor [-]

    shading_area = polyshape();
    for pv = 1:n_pv % Calculate shadow of each panel pair
        shading_area(pv) = shade_projection(PV_pairs(pv,:)', PV_width, PV_length, PV_azimuth, PV_tilt, sun_azimuth, sun_altitude, crop_width, crop_length, crop_azimuth);
    end     
  
    SF = area(union(shading_area))/(crop_width*crop_length); %shading factor of entire crop field

end

function shading_box = shade_projection(center_point, PV_width, PV_length, PV_azimuth, PV_tilt, sun_azimuth, sun_altitude, crop_width, crop_length, crop_azimuth)
    %% 2-D Panel Points 
    panel_2D = [PV_width/2, PV_length/2, 0;
                -PV_width/2, PV_length/2, 0;
                -PV_width/2, -PV_length/2, 0;
                PV_width/2, -PV_length/2,  0]; % vertices (X,Y) = (E,N)*tilted  coordinates of pv vertices wrt to local panel coordinates

    %% Calculate the 4 panel vertices with SEZ global coordinate range
    Rx_pv = [1,0,0;
             0, cosd(PV_tilt), -sind(PV_tilt);
             0, sind(PV_tilt), cosd(PV_tilt)]; % describing rotation around the x-axis given the PV_tilt angle

    Rz_pv = [-sind(PV_azimuth), -cosd(PV_azimuth), 0;
             cosd(PV_azimuth), -sind(PV_azimuth), 0;
             0, 0, 1]; % matrix describing rotation around z-axis of 90^\circ + pv_azimuth (ie angle of pv panel normal vector from due south)

    Rz_crop = [-sind(crop_azimuth), -cosd(crop_azimuth), 0;
               cosd(crop_azimuth), -sind(crop_azimuth), 0;
               0, 0, 1]; % matrix describing rotation around z-axis of 90^\circ + crop_azimuth (ie angle of crop field normal vector from due south)

    panel_global_pv = repelem(Rz_crop*center_point,1,4) + Rz_pv*Rx_pv*panel_2D';

    % Calculate projected point on cropland
    ncrop = [0,0,1]; % normal vector of crop field (assuming cropland is flat)
    sun_vector = [cosd(sun_azimuth)*cosd(sun_altitude), sind(sun_azimuth)*cosd(sun_altitude), sind(sun_altitude)]; % solar vector

    % check n * s > 0 (i.e., sun is 'in front' of crops)
    if dot(ncrop,sun_vector) < 0
        print("Warning sun is not visible. No direct radiation.")
    end
      
    PV_vertex_projected = zeros(4,3);
    for vertex =1:4
        t = -(ncrop(1)*panel_global_pv(1,vertex) + ncrop(2)*panel_global_pv(2,vertex) + ncrop(3)*panel_global_pv(3,vertex))/(ncrop(1)*sun_vector(1) + ncrop(2)*sun_vector(2) + ncrop(3)*sun_vector(3)); % intersection parameter t

        PV_vertex_projected(vertex,:) = [panel_global_pv(1,vertex)+sun_vector(1)*t, ...
                                         panel_global_pv(2,vertex)+sun_vector(2)*t, ...
                                         panel_global_pv(3,vertex)+sun_vector(3)*t]; % point where sun ray intersecting vertex touches the crop field
    end

    PV_vertex_projected_ENZ = transpose(Rz_crop)*PV_vertex_projected'; % convert back to E-N coordinate system
    PV_x = PV_vertex_projected_ENZ(1,:);
    PV_y = PV_vertex_projected_ENZ(2,:);

    crop_shape = polyshape([0,crop_width, crop_width, 0], [0,0, crop_length,crop_length]);
    shading_box = intersect(polyshape(PV_x, PV_y), crop_shape); % Only consider shading that is within the crop field

end

