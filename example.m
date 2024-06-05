% ----------------------------------------------------------------------- %
% SYLVA, https://sylva.bioaerosol.eu/
% @author: X. Shang
% ----------------------------------------------------------------------- %

%% initial lidar data and pre-defined layer information
% ###### read lidar products
% # Read files (level 2) including BSC_355, BSC_532, PDR_532, with dimension [height, time]. 
% These parameters need to have same vertical resolution (use interpolation).
% height (or altitude) is used for polle layer definition.
% start_time is used to select time window.
file = 'LidarData_KUO.nc';
BSC_355 = ncread(file, 'backscatter_coefficient_355');
BSC_532 = ncread(file, 'backscatter_coefficient_532');
PDR_532 = ncread(file, 'particle_depolarization_532');
height = ncread(file, 'altitude');
start_time = seconds(ncread(file, 'start_time')) + datetime(1970,1,1);
stop_time = seconds(ncread(file, 'stop_time')) + datetime(1970,1,1);

% # to used the near real time dust forecase, the example lidar data time
% are manually changed from 2022 to 2024.
start_time = start_time+years(4);
stop_time = stop_time+years(4);

% # define pollen layer
layer_bot = 800;    % in [m]
layer_top = 2500;   % in [m]

% # defind time window
period_begin = datetime(2024,5,23);
period_end = datetime(2024,5,27);

% # min dust concentration for dust free conditions 
threshold_dust_cnc = 5e-9; % [kg/m3], 5 ug/m3.

% # priliminary value for BAE pure pollen
BAE_pure = 0;    % assumption of BAE of pure pollen

%% ############################## main algorithm ####################### %%
% # retrieval of characteristics

% ###### remove dust period
% # dust forecase data are from fauto_SILAM_dustfore.m
flag_dust = false(1,size(PDR_532,2));
for i_time = 1:length(start_time)
    lidar_time1 = start_time(i_time);
    if lidar_time1<period_begin || lidar_time1>period_end
        continue;
    end
    try
        silam_file = dir(['silam_dust_',char(lidar_time1,'yyyy-MM-dd'),'_Kuopio.nc']).name;

        cnc_dust = ncread(silam_file, 'cnc_dust'); %[kg/m3]
        silam_cnc_dust = permute(cnc_dust(1,1,:,:),[3,4,1,2]); % [height,time]
        clear cnc_dust
        silam_height = ncread(silam_file, 'height');
        time0 = ncread(silam_file, 'time');
        t_units = ncreadatt(silam_file,"time","units");
        t0 = datetime(t_units(13:31),'InputFormat','uuuu-MM-dd hh:mm:ss');
        silam_time = t0+hours(time0);
        clear time0 t_units t0

        silam_it1 = find(silam_time<start_time(i_time),1,"last");
        silam_it2 = find(silam_time>stop_time(i_time),1,"first");
        silam_ih1 = find(silam_height<layer_bot,1,"last");
        silam_ih2 = find(silam_height>layer_top,1,"first");
        if max(silam_cnc_dust(silam_ih1:silam_ih2,silam_it1:silam_it2),[],"all")>threshold_dust_cnc
            flag_dust(i_time)=true;
        end
    catch
        warning('no SILAM dust file');
    end
end

% ###### select time window and layer height
% pollen layer
ixLayer = find(height>=layer_bot & height<layer_top); % in [m]
% analysing period
ixTime = find(start_time>=period_begin & start_time<period_end);

% ###### PDR and BAE non-linear regression fitting
[DR_pure,a1,a2,lambda1,lambda2,input_PDR,input_BAE,input_eta] = ...
    f_PDR_BAE_fitting(PDR_532,BSC_532,BSC_355,flag_dust,ixLayer,ixTime,BAE_pure);
fprintf('DR pure pollen = %.2f\n',DR_pure)

% ###### if want a final plot
x_fit = 0:0.01:0.3;
eta_fit = (a1+a2)./(x_fit+1) - a2;
bae_fit = - log(eta_fit)./log(lambda2/lambda1);

figure;
plot(input_PDR, input_BAE,'o');hold all;
plot(x_fit,bae_fit, 'r-')
plot(DR_pure,BAE_pure,'s')
xlabel(sprintf('PDR (lambda %d)',lambda1))
ylabel(sprintf('BAE (lambda %d, %d)',lambda1,lambda2))
title(sprintf('DR for pure pollen: %.2f',DR_pure))
legend('Original data','Fitted curve','Pure values')
