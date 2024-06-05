function f_SILAM_dustfore
% ----------------------------------------------------------------------- %
% SYLVA, https://sylva.bioaerosol.eu/
% @author: X. Shang
%
% Automatic download the SILAM Forecast Model for dust & plot dust products
% Websites:
% https://silam.fmi.fi/thredds/catalog.html
% ----------------------------------------------------------------------- %

date_start = char(datetime('today'),'yyyy-MM-dd');
date_end   = char(datetime('tomorrow'),'yyyy-MM-dd');

pilot_stations = {'Potenza', 'Kuopio', 'Granada', 'Hohenpeissenberg'};
for ks = 1:length(pilot_stations)
    % # step1: download .nc file from the SILAM web
    station = pilot_stations{ks};
    switch station
        case 'Potenza'
            station_lat = 40.6000;
            station_lon = 15.7200;
            coordinate_subset = 'north=40.61&west=15.72&east=15.73&south=40.60&';
        case 'Kuopio'
            station_lat = 62.7378;
            station_lon = 27.5431;
            coordinate_subset = 'north=62.74&west=27.54&east=27.55&south=62.73&';
        case 'Granada'
            station_lat = 37.164;
            station_lon = -3.605;
            coordinate_subset = 'north=37.17&west=-3.61&east=-3.60&south=37.16&';
        case 'Hohenpeissenberg'
            station_lat = 47.8019;
            station_lon = 11.0119;
            coordinate_subset = 'north=47.81&west=11.01&east=11.02&south=47.80&';
    end

    % # set the website
    silam_web = 'https://silam.fmi.fi/thredds/ncss/';
    silam_version = 'silam_europe_v5_9';
    web_subs = '/runs/';
    tmp1 = '_RUN_';
    time_start = 'T00:00:00Z?';
    % # selected variables, can be changed.
    silam_var = 'var=BLH&var=dust_tropcol&var=ocd_dust_m1_5_w550&var=ocd_dust_m20_w550&var=ocd_dust_m6_0_w550&var=ocd_dust_m_30_w550&var=ocd_dust_w550&var=air_dens&var=cnc_dust&var=cnc_dust_m1_5&var=cnc_dust_m20&var=cnc_dust_m6_0&var=cnc_dust_m_30&var=pressure&var=temperature&';
    tmp2 = 'disableProjSubset=on&horizStride=1&';
    time_window = ['time_start=',date_start,'T01%3A00%3A00Z&','time_end=',date_end,'T00%3A00%3A00Z&timeStride=1&vertCoord=&accept=netcdf'];
    url = [silam_web,silam_version,web_subs,silam_version,tmp1,date_start,time_start,silam_var,coordinate_subset,tmp2,time_window];
    % print(url)
    % this is a NCSS Request URL. copy this url and open it with a browser (e.g. chrom), the .nc file can be downloaded directly.
    
    % # download .nc file for the url
    silam_file = ['silam_dust_',date_start,'_',station,'.nc'];
    try
        websave(silam_file,url,weboptions("Timeout",50));
        fprintf('File downloaded successfully: %s\n',station)
    catch
        try
            date_start = datestr(today-1,'yyyy-mm-dd');
            url = [silam_web,silam_version,web_subs,silam_version,tmp1,date_start,time_start,silam_var,coordinate_subset,tmp2,time_window];
            websave(silam_file,url,weboptions("Timeout",50));
            fprintf('File downloaded successfully (2nd try): %s\n',station')
        catch
            warning('Failed to download file: %s\n',station)
            continue;
        end
    end

    % # step2: read and plot dust products
    time0 = ncread(silam_file, 'time');
    t_units = ncreadatt(silam_file,"time","units");
    t0 = datetime(t_units(13:31),'InputFormat','uuuu-MM-dd hh:mm:ss');
    t1_utc = t0+hours(time0);
    time = t1_utc;
    height = ncread(silam_file, 'height');
    lat = ncread(silam_file, 'lat');
    if length(lat)>1
        [~,i_lat] = min(abs(lat-station_lat));
    else
        i_lat = 1;
    end
    lat = lat(i_lat);
    lon = ncread(silam_file, 'lon');
    if length(lon)>1
        [~,i_lon] = min(abs(lon-station_lon));
    else
        i_lon = 1;
    end
    lon = lon(i_lon);
    cnc_dust = ncread(silam_file, 'cnc_dust');
    ocd_dust_w550 = ncread(silam_file, 'ocd_dust_w550');
    data_cnc_dust = permute(cnc_dust(1,1,:,:),[3,4,1,2]);
    data_ocd_dust_w550 = permute(ocd_dust_w550(1,1,:),[3,1,2]);
    data_cnc_dust_w550_1st = data_cnc_dust(1,:);
    data_cnc_dust_w550_to3km = sum(data_cnc_dust(1:7,:),1);
    data_cnc_dust_w550_all = sum(data_cnc_dust,1);

    f1 = figure;
    pt = tiledlayout(3,1,'TileSpacing','compact','Padding','loose');
    ax(1) = nexttile(1);
    plot(time,data_ocd_dust_w550,'s-');
    grid on;
    title('optical column depth dust @ 550nm');
    ylabel('Duct OD')
    ax(2) = nexttile(2);
    plot(time,data_cnc_dust_w550_all,'s-');
    grid on;hold all;
    plot(time,data_cnc_dust_w550_to3km,'o-');
    plot(time,data_cnc_dust_w550_1st,'x-');
    legend('all','0-3km','1st bin','Location','best');
    title('Concentration in air dust: SUM');
    ylabel('Layer conc.');
    ax(3) = nexttile(3);
    pcolor(time,height,data_cnc_dust);shading flat;
    ylabel('height agl [m]');
    colorbar;
    title('Concentration in air dust [kg/m^3] (0.30 + 1.5 + 6.0 + 20) \mum');
    ylim([0 3000])
    linkaxes(ax,'x');
    xlabel('UTC')
    title(pt,sprintf('%s, lat %.2f, lon %.2f (%s)',strrep(silam_version,'_','.'),lat,lon,station) )

    print('-dpng','-r300',[strrep(silam_file,'.nc', '.png')] )

    close(f1)
end