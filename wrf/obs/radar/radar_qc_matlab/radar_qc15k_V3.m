
% PARAMETROS A MODIFICAR
%*********************************************************************
%Folders

BASEDIR = '/home/juan.guerrieri/datosmate/20091117/';
SODIR   = '/home/jruiz/datos/DATA/OBS/LETKF_SO15KM_V3/';

%Reflectivity QC
nprocs = 8;           %Number of processors used to calculate echotop 3D
minref = -0.1;        %Minimum reflectivity (dBZ)
att_thld = 0.2;       %Attenuation factor threshold
z_thld = 5;           %Echo top computation threshold (dBZ)
snx = 5;              %Sizes of the local box to compute mean
sny = 5;
echotop_thld = 3500;  %Echo top threshold (m)
keep_close_ranges=false; 
%rhohv_thld = 0.7;     %RhoHV threshold

%Wind QC
nxb = 2;             %Sizes of the outer box to compute mean
nyb = 2;
nzb = 0;
dealiastr=5.0;       %Threshold for dealiasing edges detection

nx1 = 1;             %Local distance based filter.
ny1 = 1;
nz1 = 0;
tr1 = 4; 

nx2 = 10;             %Larger scale distance based filter.
ny2 = 10;
nz2 = 0;
tr2 = 15;

nfilterpass = 2;     %Number of times to apply the boxmean smoothing

nxws=2;              %Wind speckle filter
nyws=2;
nzws=0;
trws=0.2

low_wind_tr=1.5;     %Winds below this threshold will be eliminated.

%Superobbing
dx = 15000;          %Increments for cartesian grid (m)
dz = 500;           %Vertical increments for cartesian grid (m)
maxz = 18e3;        %Maximum height for the observations
maxrange=240e3;     
undef = -9999;      %Missing value
error_ref = 5;      %Reflectivity error in DBZ.
error_dv = 2;       %Doppler velocity error in m/s
id_ref_obs = 4001;  %Reflectivity id number
id_dv_obs = 4002;   %Doppler velocity id number
radar_type = 10;    %Radar wavelength
super_obbing_minn=2;  %Minimum number of data points in each super obbing grid box
super_obbing_dv_minrange=(3/2)*dx; %Doppler observations closser than this range will be ignored (recommended value (3/2)*dx)
                                   %3/2 * dx subestimation of maximum wind speeds should be less than 3%
sowindowl=900;  %Length of the superobbing window.


%*********************************************************************

initialized_so=false;

%Start matlab parallel tasks
matlabpool(nprocs)

system( ['mkdir ' SODIR] );

%List all files to read
filelist = dir( [BASEDIR '/*.nc3'] );

%Loop over files
for ifile=1:size(filelist,1)
    %Read variables in file "filename"
    filename = [ BASEDIR '/' filelist(ifile).name ];
    [radar, ref, dv, dvc, phidp, rhohv, kdp, present_ref, present_dv, present_dvc, present_phidp, present_rhohv, present_kdp] = read_radar_netcdf(filename , minref );
    

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
    compute_multiple_echotop_par(slice_radar,ref(1,:,:),z_thld,[],[]);

    parfor ia=1:radar.na
      [echo_top_3d(ia,:,:), ~, echo_depth_3d(ia,:,:), ~, max_dbz_z_3d(ia,:,:)] = ...
      compute_multiple_echotop_par(slice_radar, ref(ia,:,:), z_thld, televation, trange);
    end
    
    %%Perform box averaging smoothing
    [echo_top_3d] = compute_boxmax(echo_top_3d, snx, sny, 0);

    echo_top_3d( isnan(echo_top_3d) ) = 0.0;

    %Busco la altura del maximo angulo de elevacion
    tmpz=radar.Z(1,:,end);
    minrange=1;
     for ii=1:length(tmpz)-1
        if( tmpz(ii) < echotop_thld & tmpz(ii+1) > echotop_thld )
          minrange = ii;
        end
     end
    %Protect the pixels close to the radar from echo top filter
    if( keep_close_ranges )

      %Busco la altura del maximo angulo de elevacion
      tmpz=radar.Z(1,:,end);
      minrange=1;
      for ii=1:length(tmpz)-1
         if( tmpz(ii) < echotop_thld & tmpz(ii+1) > echotop_thld )
           minrange = ii;
         end
      end

       echo_top_3d(:,1:minrange,:)=echotop_thld+1.0;
    else
       %We remove all the pixels where echo top is too low.
       refc(:,1:minrange,:)=NaN;
    end

    %%Remove gates with low echo top
    refc( echo_top_3d < echotop_thld & ref > 0 ) = NaN;

    %Remove strongly attenuated beams
    [attenuation_factor] = compute_atenuation_qc(radar, ref);
    refc( attenuation_factor < att_thld ) = NaN;

    %Remove gates associated with small doppler velocity
    if ( present_dv )
       refc( abs(dv) < low_wind_tr ) = NaN;
    end

    %Remove pixels with reflectivity for elevations 1 or 2 but not over that.
    %This removes clutter generated by distant anomalous propagation.

    tmpref=refc(:,:,2);

    tmpref( tmpref > minref & ( refc(:,:,3) == minref | isnan(refc(:,:,2) ) ) )=NaN;

    refc(:,:,2)=tmpref;

    tmpref=refc(:,:,1);

    tmpref( tmpref > minref & ( refc(:,:,2) == minref | isnan(refc(:,:,2) )  ) )=NaN;

    refc(:,:,1)=tmpref;


    %Remove echoes associated with low RhoHV
    %if ( present_rhohv )
    %   rhohv( rhohv < 0 ) = NaN;
    %   [rhohvsmooth] = compute_boxmean(rhohv, 1, 1, 1);
    %   refc( rhohvsmooth < rhohv_thld ) = NaN;
    %end

    %Detect missing values within clouds
    [refc] = compute_detectmissing(radar,refc,10.0,minref,15);

    %Use texture
    %[texture]= compute_texture(radar,refc,0,3,0);
    %refc( texture > 100.0 )=NaN;

    reftmp=refc;
 
    npass=6;
       
    for ip=1:npass

      [distref]=compute_distance(refc,reftmp,1,1,0);
        
      reftmp=refc;
  
      reftmp( distref > 15.0 )=NaN;
      
     end  
    
     refc( (distref > 15.0 | isnan(distref)) & refc < 30.0 )=NaN;
%        
% 
     [count]=compute_boxcount(refc,2,2,0,minref);
     
     refc(count < 0.2 & refc > minref )=NaN;
       
    
     %Remove gates with low reflectivity below 3Km (PBL echoes)
     refc( (refc > minref & refc < 15.0) & radar.Z < 2.5e3 )=NaN;


    % WIND QC
    %*****************************************************************
    %Dealiasing is perfomed using pyart routines (region_dealiased)
    %Some faulty points may remain which are cleaned with the following functions
    if ( present_dvc )

       %Remove noise in the doppler velocity
       dvc( abs(dvc) < low_wind_tr ) = NaN;

       [ dvcc ] = compute_detectdealiasingborders( radar , dv , dvc , nxb , nyb , nzb , dealiastr );

       %npass=2;
    
       dvtmp=dvcc;
    
       for ip=1:nfilterpass

         [distdv]=compute_distance(dvcc,dvtmp,nx1,ny1,nz1);

         dvtmp=dvcc;
    
         dvtmp( distdv > tr1 )=NaN;

       end
    
       dvcc(distdv > tr1 | isnan(distdv) )=NaN;
      
       dvtmp=dvcc;

    
       for ip=1:nfilterpass

        [distdv]=compute_distance(dvcc,dvtmp,nx2,ny2,nz2);

        dvtmp=dvcc;
    
        dvtmp(distdv > tr2 )=NaN;
      
       end  
    
       dvcc(distdv > tr2 | isnan(distdv) )=NaN;

       [count]=compute_boxcount(abs(dvcc),nxws,nxws,nxws,low_wind_tr);
    
       dvcc(count < trws )=NaN;
    
% 
     else
        dvcc = NaN( size(dvc) );
     end

    %Save QC variables in radar geometry to .mat to plot in python
    latitude = radar.latitude;
    longitude = radar.longitude;
    height = radar.height;
    %file2write = [QCDIR '/qc_' filelist(ifile).name(7:21) '.mat'];
    %save(file2write, 'dvcc', 'refc', 'dv', 'dvc_ini', 'ref_ini', 'latitude', 'longitude', 'height');
    disp( '***********************************************************' );
    disp( ['DONE RADAR QC for file ' filelist(ifile).name(7:21)] );
    disp( '***********************************************************' );
    %SUPEROBBING
    %*****************************************************************
    %Generate cartesian grid
    [cart] = define_cartesian_grid_v2(radar, dx, dz, maxz , maxrange ) ;

    nx=size(cart.z,1);
    ny=size(cart.z,2);
    nz=size(cart.z,3);

    disp( ['The size of the super obbing domain is nx=' num2str(nx) ' ny=' num2str(ny) ' nz=' num2str(nz) ] );

    %Obtain date and time and add them to the radar structure 
    date_rad_ini_str = radar.time_coverage_start([1:4  6:7 9:10 12:13 15:16 18:19]);
    date_rad_ini_num = datenum( date_rad_ini_str , 'yyyymmddHHMMSS');

    %Obtain date and time and add them to the radar structure 
    date_rad_end_str = radar.time_coverage_end([1:4  6:7 9:10 12:13 15:16 18:19]);
    date_rad_end_num = datenum( date_rad_end_str , 'yyyymmddHHMMSS');

    date_rad_mean_num= 0.5* date_rad_ini_num + 0.5*date_rad_end_num  ;
    date_rad_mean_str= datestr( date_rad_mean_num ,  'yyyymmddHHMMSS' );

    display(['Current volume starts at:' date_rad_ini_str ])
    display(['Current volume ends   at:' date_rad_end_str ])
    

    if ( ~ initialized_so ) ;

       dateref=[date_rad_ini_str(1:4) '0101000000' ];
       %This is the first volume. 
       %Get the date of the first SO.
       refdate=round( date_rad_ini_num )+1.0; %Uso una fecha de referencia para evitar problemas de redondeo.
       tmp = date_rad_ini_num - refdate;
       %Fechas relativas a la refdate.
       date_so_mean_num=round(tmp*86400.0/sowindowl)*sowindowl/86400.0;
       date_so_ini_num=tmp-sowindowl/(2*86400);
       date_so_end_num=tmp+sowindowl/(2*86400);

       initialized_so = true ;

       display(['Starting a new SO window :' datestr( round( (date_so_mean_num + refdate)*86400)/86400 ,'yyyy-mm-dd-HH:MM:SS') ])

       %Esta funcion solo llena de ceros la estructura gridrad
       gridrad=initialize_gridrad_fun(cart);

    end 

    cont=true;

    while cont 

      tmpref=refc;
      tmpdv =dvcc;

      so_ini_time = 86400 * ( date_so_ini_num + refdate - date_rad_ini_num ); %Superobbing bin initial time in seconds with respect to superobbing bin start.
      so_end_time = 86400 * ( date_so_end_num + refdate - date_rad_ini_num ); %Superobbing bin end time in seconds with respect to superobbing bin start.
        
      for ii=1:radar.na
        for jj=1:radar.ne
            if( ~ (radar.time(ii,jj) <= so_end_time & radar.time(ii,jj) >= so_ini_time ) )
              %Do not include this data in the current so bin.
              tmpref(ii,:,jj) = NaN;
              tmpdv(ii,:,jj)  = NaN;
            end
        end
      end

      display(['A total number of ' num2str( (sum(sum(~isnan(tmpref))) )) ' valid reflectivity observations available from this volume to this slot'])
      display(['A total number of ' num2str( (sum(sum(~isnan(tmpdv))) )) ' valid Doppler velocity observations available from this volume to this slot'])

      [gridrad] = radar_superobbing_4d( radar , tmpref , tmpdv , cart , gridrad );

      if(   date_rad_end_num < date_so_end_num + refdate )

         cont = false ;
         %Lets keep the same so bin and get a new radar data file.

      else

         %We need to start a new superobbing window and write the previous one.
         display(['Writing data corresponding to superobbing slot centered at:' datestr( round( (date_so_mean_num + refdate)*86400)/86400 , 'yyyy-mm-dd-HH:MM:SS') ])
         fileout = [SODIR '/radar01_' datestr( round( (date_so_mean_num + refdate)*86400)/86400 ,'yyyymmddHHMMSS') ];

         gridrad.grid_dv( gridrad.grid_count_dv < super_obbing_minn )=0.0;
         gridrad.grid_count_dv( gridrad.grid_count_dv < super_obbing_minn )=0.0;
         gridrad.grid_ref( gridrad.grid_count_ref < super_obbing_minn )=0.0;
         gridrad.grid_count_ref( gridrad.grid_count_ref < super_obbing_minn )=0.0;

         %Remove super obbing observations which are close to the radar. Because azimuth chagens too much 
         %within this ranges.
         %TODO esto se puede solucionar haciendo un local VAD fit y estimando el Vr que se deberia medir
         %en el azimuth medio y la elevacion media de acuerdo con un ajuste de los datos que caen dentro del volumen.
         gridrad.grid_count_dv( gridrad.grid_ra_dv < super_obbing_dv_minrange )=0.0;
         gridrad.grid_dv( gridrad.grid_ra_dv < super_obbing_dv_minrange )=0.0;

         %Write the current so bin.
         write_so( radar , gridrad , fileout , true ); 

         %Define the next so bin.
         date_so_mean_num=date_so_mean_num + sowindowl/86400;
         date_so_ini_num=date_so_ini_num + sowindowl/86400;
         date_so_end_num=date_so_end_num + sowindowl/86400;

         display(['Starting a new SO window :' datestr( round( (date_so_mean_num + refdate)*86400)/86400 ,'yyyy-mm-dd-HH:MM:SS') ])

         %Esta funcion solo llena de ceros la estructura gridrad
         gridrad=initialize_gridrad_fun(cart);
 
      end

    end

    disp( '***********************************************************' );
    disp( ['DONE SUPEROBBING for file ' filelist(ifile).name(7:21)] );
    disp( '***********************************************************' );

end  %End loop over files

%End matlab parallel tasks
matlabpool close

disp('DONE QC AND SUPPEROBBING');


