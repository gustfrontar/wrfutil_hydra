clear all
close all

% PARAMETROS A MODIFICAR
%*********************************************************************
%Folders
OBS = 'ANGUIL';                  %Radar folder
CASE = '20100111';               %Case study folder
DATA = '240';                    %Range folder
INPUT = 'region_dealiased_data'; %Radar data folder
QC = 'QC_DATA';                  %QC folder
SO = 'LETKF_SO2KM';              %Supperobing folder

%Radar
undef = -9999;      %Missing value
error_ref = 5;      %Reflectivity error in DBZ.
error_dv = 2;       %Doppler velocity error in m/s
id_ref_obs = 4001;  %Reflectivity id number
id_dv_obs = 4002;   %Doppler velocity id number
radar_type = 10;    %Radar wavelength

%Reflectivity QC
nprocs = 8;           %Number of processors used to calculate echotop 3D
minref = 0;           %Minimum reflectivity (dBZ)
att_thld = 0.2;       %Attenuation factor threshold
z_thld = 5;           %Echo top computation threshold (dBZ)
snx = 5;              %Sizes of the local box to compute mean
sny = 5;
echotop_thld = 3000;  %Echo top threshold (m)
rhohv_thld = 0.8;     %RhoHV threshold

%Wind QC
npixels = 3;         %Number of pixels to remove next to dealiased region's edges
nx = 10;             %Sizes of the outer box to compute mean
ny = 10;
nz = 0;
nx2 = 4;             %Sizes of the inner box to compute mean
ny2 = 4;
nz2 = 0;
vad_thld = 25;       %VAD computation threshold
nanomaly_thld = 3.0; %Normalized anomaly velocity threshold
nfilterpass = 1;     %Number of times to apply the boxmean smoothing

%Superobbing
dx = 2000;    %Increments for cartesian grid (m)
dz = 2000;
maxz = 25e3;  %Maximum height for the observations

%*********************************************************************

%Start matlab parallel tasks
matlabpool(nprocs)

%Create necessary directories
BASEDIR = '/home/paula.maldonado/datosmate/RADAROBS';
FILEDIR = [BASEDIR '/' OBS '/' CASE '/' DATA '/' INPUT];
QCDIR = [FILEDIR '/' QC];
SODIR = [FILEDIR '/' SO];
system( ['mkdir ' QCDIR] );
system( ['mkdir ' SODIR] );

%List all files to read
filelist = dir( [FILEDIR '/*.nc3'] );

%Loop over files
for ifile=1:size(filelist,1)
    %Read variables in file "filename"
    filename = [ FILEDIR '/' filelist(ifile).name ];
    [radar, ref, dv, dvc, phidp, rhohv, kdp, present_ref, present_dv, present_dvc, present_phidp, present_rhohv, present_kdp] = read_radar_netcdf(filename);
 
    %Save reflectivity and doppler velocity with no QC (obs: dvc is dealiased)
    ref_ini = ref;
    dvc_ini = dvc;

    %Georeference radar volume. It will be performed for each time due to possible changes
    %in observation strategy
    [radar] = georeference_radar_data(radar);

    %Add features to radar structure
    radar.replacerefmissing = minref;
    radar.error_ref = error_ref;
    radar.error_dv = error_dv;
    radar.id_ref_obs = id_ref_obs;
    radar.id_dv_obs = id_dv_obs;
    radar.radar_type = radar_type;
    radar.undef = undef;

    %Topography blocking.
    % TODO: En la version de python-fortran debe estar incorporado este
    % chequeo (importante sobre todo para radares como el de Cordoba).

    % REFLECTIVITY QC
    %*****************************************************************
    %Set reflectivity minimum
    ref( ref < minref ) = minref;
    refc=ref;

    %Remove strongly attenuated beams
    [attenuation_factor] = compute_atenuation_qc(radar, ref);
    refc( attenuation_factor < att_thld ) = NaN;

    %Remove gates associated with small doppler velocity
    if ( present_dv )
       refc( abs(dv) < 1 ) = NaN;
    end

    %Remove gates with low echo top (clutter, PBL echoes)
    %%Define necessary varibles
    echo_top_3d = NaN(radar.na, radar.nr, radar.ne);
    echo_depth_3d = NaN(radar.na, radar.nr, radar.ne);
    max_dbz_z_3d = NaN(radar.na, radar.nr, radar.ne);
    echo_max_dbz_z = NaN(radar.na, radar.nr, radar.ne);

    slice_radar.na = 1;
    slice_radar.nr = radar.nr;
    slice_radar.Rh = radar.Rh;
    slice_radar.ne = radar.ne;
    slice_radar.elevation = radar.elevation;
    slice_radar.height = radar.height;
    slice_radar.replacerefmissing = radar.replacerefmissing;
    slice_radar.range = radar.range;
    slice_radar.Z = radar.Z(1,:,:);
    slice_radar.beam_wid_v = radar.beam_wid_v;
    slice_radar.beam_wid_h = radar.beam_wid_h;

    %%Compute echo_top with a local definition (each cloud has its own
    %%echo top and we can have several echo tops at the same horizontal
    %%location correspoding for overlaping clouds)
    [~, ~, ~, ~, ~, televation, trange] = ...
    compute_multiple_echotop_par(slice_radar, ref(1,:,:), z_thld, [], []);

    parfor ia=1:radar.na
      [echo_top_3d(ia,:,:), ~, echo_depth_3d(ia,:,:), ~, max_dbz_z_3d(ia,:,:)] = ...
      compute_multiple_echotop_par(slice_radar, ref(ia,:,:), z_thld, televation, trange);
    end
    
    %%Perform box averaging smoothing
    [echo_top_3d] = compute_boxmean(echo_top_3d, snx, sny, 0);
    [echo_depth_3d] = compute_boxmean(echo_depth_3d, snx, sny, 0);
    [max_dbz_z_3d] = compute_boxmean(max_dbz_z_3d, snx, sny, 0);

    echo_top_3d( isnan(echo_top_3d) ) = 0.0;
    echo_depth_3d( isnan(echo_depth_3d) ) = 0.0;
    max_dbz_z_3d( isnan(max_dbz_z_3d) ) = 0.0;

    %%Remove gates with low echo top
    refc( echo_top_3d < echotop_thld & ref > 0 ) = NaN;

    %Remove echoes associated with low RhoHV
    if ( present_rhohv )
       rhohv( rhohv < 0 ) = NaN;
       [rhohvsmooth] = compute_boxmean(rhohv, 1, 1, 1);
       refc( rhohvsmooth < rhohv_thld ) = NaN;
    end

    %Remove gates with low reflectivity below 3Km (PBL echoes)
    %refc( refc < 25 & radar.Z < 3e3 )=NaN;

    % WIND QC
    %*****************************************************************
    %Dealiasing is perfomed using pyart routines (region_dealiased)
    %Some faulty points may remain which are cleaned with the following functions
    if ( present_dvc )

       %Remove noise in the doppler velocity
       dvc( abs(dvc) < 1 ) = NaN;

       %Remove dealiased regions edges
       diff = dv-dvc_ini;
       
       %%Loop over elevations
       for k=1:radar.ne

         %%Loop over azimuths
         for i=npixels:radar.na-npixels

           %%Loop over ranges
           for j=npixels:radar.nr-npixels
             if ( diff(i,j,k) ~= 0 && ~isnan(diff(i,j,k)) )
                pixel = diff(i,j,k);

                %%Loop over pixel's next neighbours
                ineighbour = -1;
                while ( ineighbour < 2 )

                   %%Loop over pixels to be remove
                   ipixels = 1;
                   while ( ipixels <= npixels )

                      neighbour_az = diff(i+ineighbour,j,k);
                      if ( abs(pixel-neighbour_az) > 70 )
                         dvc(i+ineighbour*ipixels,j,k) = NaN;
                      end
                      neighbour_ra = diff(i,j+ineighbour,k);
                      if ( abs(pixel-neighbour_ra) > 70 )
                         dvc(i,j+ineighbour*ipixels,k) = NaN;
                      end
                      ipixels = ipixels + 1;

                   end  %End loop over pixels to remove
                   ineighbour = ineighbour + 2;

                end  %End loop over pixel's next neighbour

             end  
           end  %End loop over ranges

         end  %End loop over azimuths

       end  %End loop over elevations

       %Remove outliers using VAD routine
       [dvcc, u_vad, v_vad, uo, vo, tmp_wind, vr_background, sample_coverage] = compute_wind_qc(radar, dvc, nx, ny, nz, nx2, ny2, nz2, vad_thld, nanomaly_thld);

       %Remove outliers using box avreaging smoothing with different box sizes
       for ifilter=1:nfilterpass
         [mean_dvcc, std_dvcc] = compute_boxmean_V2(dvcc, nx, ny, nz, nx2, ny2, nz2);
         std_dvcc ( std_dvcc == 0 ) = NaN;
         nanomaly_dvcc = abs(dvcc-mean_dvcc)./std_dvcc;
         nanomaly_dvcc ( nanomaly_dvcc > nanomaly_thld ) = NaN; 
         dvcc ( isnan(nanomaly_dvcc) ) = NaN;
       end
      
       for ifilter=1:nfilterpass
         [mean_dvcc, std_dvcc] = compute_boxmean_V2(dvcc, round(nx/2), round(ny/2), round(nz/2), round(nx2/2), round(ny2/2), round(nz2/2));
         std_dvcc ( std_dvcc == 0 ) = NaN;
         nanomaly_dvcc = abs(dvcc-mean_dvcc)./std_dvcc;
         nanomaly_dvcc ( nanomaly_dvcc > 0.8*(nanomaly_thld) ) = NaN;
         dvcc ( isnan(nanomaly_dvcc) ) = NaN;
       end

    else
       dvcc = NaN( size(dvc) );
    end

    %Save QC variables in radar geometry to .mat to plot in python
    latitude = radar.latitude;
    longitude = radar.longitude;
    height = radar.height;
    file2write = [QCDIR '/qc_' filelist(ifile).name(7:21) '.mat'];
    save(file2write, 'dvcc', 'refc', 'dv', 'dvc_ini', 'ref_ini', 'latitude', 'longitude', 'height');
    disp( '***********************************************************' );
    disp( ['DONE RADAR QC for file ' filelist(ifile).name(7:21)] );
    disp( '***********************************************************' );
    %SUPEROBBING
    %*****************************************************************
    %Generate cartesian grid
    [cart] = define_cartesian_grid(radar, dx, dz, maxz) ;

    %Obtain date and time and add them to the radar structure 
    date = radar.time_coverage_start([1:4  6:7 9:10]);
    hour = str2num(radar.time_coverage_start(12:13));
    min = str2num(radar.time_coverage_start(15:16));
    sec = str2num(radar.time_coverage_start(18:19));

    radar.year = str2num(date(1:4));
    radar.month = str2num(date(5:6));
    radar.day = str2num(date(7:8));
    radar.hour = hour;
    radar.minute = min;
    radar.second = sec;

    %Set other parameters in radar structure
    radar.beam_wid_h = 1;
    radar.beam_wid_v = 1;
    radar.attenuation_factor = 0;  %Global attenuation factor (i.e. rain over the radome)

    %Loop over elevations
    minelev = 1;
    for ii=1:radar.ne  
       if ( ii > 1 )
          %Define strings to create the fileout name
          strhour = num2str(hour);
          if ( hour < 10 ); strhour = ['0' strhour]; end
          strmin = num2str(min);
          if ( min < 10 ); strmin = ['0' strmin]; end
          %strsec = num2str(sec);
          %if ( sec < 10 ); strsec = ['0' strsec]; end
          fileout = [SODIR '/radar01_' date strhour strmin '00.dat'];
          %Use the time corresponding to the first azimuth to decide if
          %the scan is going to be included or not in the current group
          sec = sec + ( radar.time(1,ii) - radar.time(1,ii-1) );
       end

       %Divide the radar volume into 1min files 
       if ( sec > 60 )
          %fileout
          %Determine last elevation to include in the corresponding file
          endelev = ii - 1;
          
          %Define a temporal radar structure
          tmp_radar = radar;
          tmp_radar.elevation = radar.elevation(minelev:endelev);
          tmp_radar.ne = ( endelev - minelev ) + 1;
          tmp_radar.Z = radar.Z(:,:,minelev:endelev);
          tmp_radar.longitude = radar.longitude(:,:,minelev:endelev);
          tmp_radar.latitude = radar.latitude(:,:,minelev:endelev);
          tmp_radar.Rh = radar.Rh(:,minelev:endelev);
          tmp_radar.height = radar.height(:,minelev:endelev);

          %Perform superobbing
          [grid_ref, grid_count_ref, grid_dv, grid_count_dv, grid_az_ref, grid_el_ref, grid_ra_ref] = ...
          radar_superobbing(tmp_radar, refc(:,:,minelev:endelev), dvcc(:,:,minelev:endelev), cart, fileout);

          minelev = endelev + 1; 

          sec = sec - 60;
          min = min + 1;

          if ( min > 60 )
             hour = hour + 1;
             min = min - 60;
          end

          if ( hour > 23 )
             hour = 0;
             date = datestr( datenum(date, 'yyyymmdd') + 1, 'yyyymmdd' );
          end
       end

    end  %End over elevations

    %Write the results in a binary data structure (For analysis and
    %forecast verification).
    fileout = [SODIR '/radar01_' date strhour strmin '00.rad'];
    qcflag = zeros(size(ref_ini));
    write_radar_data_seq(radar, refc, dvcc, attenuation_factor, qcflag, fileout, 'b');
    disp( '***********************************************************' );
    disp( ['DONE SUPEROBBING for file ' filelist(ifile).name(7:21)] );
    disp( '***********************************************************' );

end  %End loop over files

%End matlab parallel tasks
matlabpool close

disp('DONE QC AND SUPPEROBBING');


%EXTRA: to plot qc and so results
%********************************************************************
%figure
%subplot(2,2,1)
%  pcolor(longitude(:,:,1),latitude(:,:,1),ref_ini(:,:,8))
%  shading flat;
%  caxis([-20 70])
%subplot(2,2,2)
%  pcolor(longitude(:,:,1),latitude(:,:,1),refc(:,:,8))
%  shading flat;
%  caxis([-20 70])
%subplot(2,2,3)
%  pcolor(longitude(:,:,1),latitude(:,:,1),dvc_ini(:,:,8))
%  shading flat;
%  caxis([-50 50])
%subplot(2,2,4)
%  pcolor(longitude(:,:,1),latitude(:,:,1),dvcc(:,:,8))
%  shading flat;
%  caxis([-50 50])


%ref_dbz=10.*log10(grid_ref);
%figure
%subplot(4,2,1)
%  pcolor(longitude(:,:,3),latitude(:,:,3),refc(:,:,2));
%  shading flat;colorbar
%  caxis([0 70])
%subplot(4,2,2)
%  pcolor(longitude(:,:,3),latitude(:,:,3),dvcc(:,:,2));
%  shading flat;colorbar
%  caxis([-50 50])
%subplot(4,2,3)
%  pcolor(cart.lon,cart.lat,ref_dbz(:,:,2));shading flat;colorbar
%  caxis([0 70])
%subplot(4,2,4)
%  pcolor(cart.lon,cart.lat,grid_dv(:,:,2));shading flat;colorbar
%  caxis([-50 50])
%subplot(4,2,5)
%  pcolor(cart.lon,cart.lat,ref_dbz(:,:,6));shading flat;colorbar
%  caxis([0 70])
%subplot(4,2,6)
%  pcolor(cart.lon,cart.lat,grid_dv(:,:,6));shading flat;colorbar
%  caxis([-50 50])
%subplot(4,2,7)
%  pcolor(cart.lon,cart.lat,ref_dbz(:,:,13));shading flat;colorbar
%  caxis([0 70])
%subplot(4,2,8)
%  pcolor(cart.lon,cart.lat,grid_dv(:,:,13));shading flat;colorbar
%  caxis([-50 50])


