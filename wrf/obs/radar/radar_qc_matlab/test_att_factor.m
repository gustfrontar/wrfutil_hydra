clear all
close all

% PARAMETROS
OBS = 'ANGUIL';
CASE = '20100111';
DATA = '120';
INPUT = 'region_dealiased_data';
nprocs=8; %Numero de procesadores que vamos a usar para calcular echo top 3D.
%************************************************

BASEDIR = '/home/paula.maldonado/datosmate/RADAROBS';
FILEDIR = [ BASEDIR '/' OBS '/' CASE '/' DATA '/' INPUT ];
filelist=dir([FILEDIR '/*.nc3']);

undef=-9999;
snx=5;
sny=5;
error_ref=5; %Reflectivity error in DBZ.
error_dv=2;  %Doppler velocity error in m/s.
id_ref_obs=4001;
id_dv_obs=4002;
radar_type=1;

%Increments for cartesian grid (m).
dx=2000;
dz=2000;
maxz=25e3;  %Maximum height for the observations.

%QC
minref=0;
zthresh=5; %Threshold in dBZ for echo top computation.
rhohvt=0.8; %Rhohv threshold for QC
echo_top_t=3000;%Echo top threshold (m).
attenuation_t=0.4;%Attenuation factor threshold.


for ifile=1:size(filelist,1)

    filename=[ FILEDIR '/' filelist(ifile).name ]
    [radar , ref , dv , dvc , phidp , rhohv , kdp , present_ref , present_dv , present_dvc , present_phidp , present_rhohv, present_kdp ]=read_radar_netcdf(filename);

    ref_ini=ref;
    dvc_ini=dvc;
    qcflag=zeros(size(ref));

    %Georeference will be performed for each time due to possible changes
    %in observation strategy.
    [radar]=georeference_radar_data(radar);
    radar.replacerefmissing=minref;
    radar.error_ref = error_ref;
    radar.error_dv = error_dv ;
    radar.id_ref_obs = id_ref_obs;
    radar.id_dv_obs = id_dv_obs;
    radar.radar_type = radar_type;
    radar.undef = undef;

    %Topography blocking.
    % TODO: En la version de python-fortran debe estar incorporado este
    % chequeo (importante sobre todo para radares como el de Cordoba).

    % REFLECTIVITY QC

    ref(ref < minref)=minref;
    refc=ref;
    refc_01=ref;
    refc_02=ref;
    refc_03=ref;
    refc_04=ref;
    refc_06=ref;


    %Remove strongly attenuated beams
    [attenuation_factor]=compute_atenuation_qc(radar,ref);
    refc_01( attenuation_factor < 0.1 )=NaN;
    refc_02( attenuation_factor < 0.2 )=NaN;
    refc_03( attenuation_factor < 0.3 )=NaN;
    refc_04( attenuation_factor < 0.4 )=NaN;
    refc_06( attenuation_factor < 0.6 )=NaN;

    %Remove gates associated with small doppler velocity.
    if( present_dv )
       refc(abs(dv) < 1)=NaN;
       refc_01(abs(dv) < 1)=NaN;
       refc_02(abs(dv) < 1)=NaN;
       refc_03(abs(dv) < 1)=NaN;
       refc_04(abs(dv) < 1)=NaN;
       refc_06(abs(dv) < 1)=NaN;
    end

    %Remove gates with low echo top (clutter, PBL echoes)

    echo_top_3d=NaN(radar.na,radar.nr,radar.ne);
    echo_depth_3d=NaN(radar.na,radar.nr,radar.ne);
    max_dbz_zi_3d=NaN(radar.na,radar.nr,radar.ne);
    %echo_max_dbz_z=NaN(radar.na,radar.nr,radar.ne);

    %Compute echo_top V2 0 with a local definition (each cloud has its own
    %echo top and we can have several echo tops at the same horizontal
    %location correspoding to overlaping clouds.

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Calculo echo top 3D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    slice_radar.na=1;
    slice_radar.nr=radar.nr;
    slice_radar.Rh=radar.Rh;
    slice_radar.ne=radar.ne;
    slice_radar.elevation=radar.elevation;
    slice_radar.height=radar.height;
    slice_radar.replacerefmissing=radar.replacerefmissing;
    slice_radar.range=radar.range;
    slice_radar.Z=radar.Z(1,:,:);
    slice_radar.beam_wid_v=radar.beam_wid_v;
    slice_radar.beam_wid_h=radar.beam_wid_h;
    %
    [~, ~ , ~ , ~ , ~ , televation , trange ]=compute_multiple_echotop_par(slice_radar,ref(1,:,:),zthresh,[],[]);

    %matlabpool(nprocs)

    %parfor ia=1:radar.na

    %[echo_top_3d(ia,:,:) , ~ , echo_depth_3d(ia,:,:) , ~ , max_dbz_zi_3d(ia,:,:) ]=compute_multiple_echotop_par(slice_radar,ref(ia,:,:),zthresh,televation,trange);

    %end

    %[echo_top_3d]=compute_boxmean(echo_top_3d,snx,sny,0);
    %[echo_depth_3d]=compute_boxmean(echo_depth_3d,snx,sny,0);
    %[max_dbz_z_3d]=compute_boxmean(max_dbz_z_3d,snx,sny,0);

    %Replace NaN, in all cases by values that increases clutter
    %probability.

    %echo_top_3d(isnan(echo_top_3d))=0.0;
    %echo_depth_3d(isnan(echo_depth_3d))=0.0;
    %max_dbz_zi_3d(isnan(max_dbz_zi_3d))=0.0;


%    [echo_top_3d ,echo_depth_3d ,echo_depth_size_3d ,total_max_dbz_3d ,total_max_dbz_level_3d ,max_dbz_3d ,max_dbz_level_3d]= ...
%    compute_echotop3d(radar,ref,zthresh,snx,sny);
    %refc(echo_top_3d < echo_top_t & ref > 0 )=NaN;
    %refc_01(echo_top_3d < echo_top_t & ref > 0 )=NaN;
    %refc_02(echo_top_3d < echo_top_t & ref > 0 )=NaN;
    %refc_03(echo_top_3d < echo_top_t & ref > 0 )=NaN;
    %refc_04(echo_top_3d < echo_top_t & ref > 0 )=NaN;
    %refc_06(echo_top_3d < echo_top_t & ref > 0 )=NaN;

    %Remove echoes associated with low RhoHV
    if( present_rhohv )
       rhohv(rhohv < 0)=NaN;
       [rhohvsmooth]=compute_boxmean(rhohv,1,1,1);
       refc( rhohvsmooth < rhohvt ) = NaN;
       refc_01( rhohvsmooth < rhohvt ) = NaN;
       refc_02( rhohvsmooth < rhohvt ) = NaN;
       refc_03( rhohvsmooth < rhohvt ) = NaN;
       refc_04( rhohvsmooth < rhohvt ) = NaN;
       refc_06( rhohvsmooth < rhohvt ) = NaN;
    end

    %Remove gates with low reflectivity below 3Km (PBL echoes)
    %refc( refc < 25 & radar.Z < 3e3 )=NaN;

    % WIND QC
    %(Dealiasing is perfomed using pyart routines, but some faulty
    %points may remain which are cleaned with the following function).

    %Remove noise in the doppler velocity
    %if( present_dvc )
       %dvc(dvc < -200) =NaN;
       %dvc(dvc == 0)= NaN;

       %Utilizar la rutina basada en el VAD para remover los puntos erroneos.
       %[dvcc , u_vad , v_vad , uo , vo , tmp_wind , vr_background , sample_coverage]=compute_wind_qc(radar,dvc,5,5,1,25,10);

       %dvccc=dvcc;

       %Remove low doppler velocity values (possible clutter).
       %dvccc(abs(dvccc) < 1)=NaN;

    %end

    latitude=radar.latitude;
    longitude=radar.longitude;

    ilev=1;
    while (ilev < 10)

    hFig=figure;
    subplot(3,2,1)
    pcolor(longitude(:,:,ilev),latitude(:,:,ilev),ref_ini(:,:,ilev))
    shading flat; colorbar
    title('RAW REF')
    caxis([-20 70])
    subplot(3,2,2)
    pcolor(longitude(:,:,ilev),latitude(:,:,ilev),refc_01(:,:,ilev))
    shading flat; colorbar
    title('ATT FACTOR 0.1')
    caxis([-20 70])
    subplot(3,2,3)
    pcolor(longitude(:,:,ilev),latitude(:,:,ilev),refc_02(:,:,ilev))
    shading flat; colorbar
    title('ATT FACTOR 0.2')
    caxis([-20 70])
    subplot(3,2,4)
    pcolor(longitude(:,:,ilev),latitude(:,:,ilev),refc_03(:,:,ilev))
    shading flat; colorbar
    title('ATT FACTOR 0.3')
    caxis([-20 70])
    subplot(3,2,5)
    pcolor(longitude(:,:,ilev),latitude(:,:,ilev),refc_04(:,:,ilev))
    shading flat; colorbar
    title('ATT FACTOR 0.4')
    caxis([-20 70])
    subplot(3,2,6)
    pcolor(longitude(:,:,ilev),latitude(:,:,ilev),refc_06(:,:,ilev))
    shading flat; colorbar
    title('ATT FACTOR 0.6')
    caxis([-20 70])

    print(['/home/paula.maldonado/datosmate/RADAROBS/ANGUIL/20100111/120/region_dealiased_data/test_att_factor/elev_' num2str(ilev) '_' filelist(ifile).name(7:21) '.png'], '-dpng')

    close(hFig)

    ilev=ilev+2;

    end


    %latitude=radar.latitude;
    %longitude=radar.longitude;
    %file2write=[ FILEDIR '/qc_radar_' filelist(ifile).name(7:21) '.mat' ];
    %save(file2write, 'dvccc', 'refc', 'dvc_ini', 'ref_ini', 'latitude', 'longitude')

    %nccreate(file2write,'Vda_corrected');
    %ncwrite(file2write,'Vda_corrected', dvccc);
    %nccreate(file2write,'dBZ_corrected');
    %ncwrite(file2write,'dBZ_corrected', refc);


    %SUPEROBBING

    %%Dividir los archivos en conjuntos de 1 minuto.
    %[ cart ] = define_cartesian_grid( radar , dx , dz , maxz) ;

    %%[tmppath tmpfilename]=fileparts(filename);

    %date=radar.time_coverage_start([1:4  6:7 9:10]);
    %hour=str2num(radar.time_coverage_start(12:13));
    %min=str2num(radar.time_coverage_start(15:16));
    %sec=str2num(radar.time_coverage_start(18:19));

    %radar.year=str2num(date(1:4));
    %radar.month=str2num(date(5:6));
    %radar.day=str2num(date(7:8));
    %radar.hour=hour;
    %radar.minute=min;
    %radar.second=sec;


    %radar.beam_wid_h=1;
    %radar.beam_wid_v=1;

    %radar.attenuation_factor=0;  %Global attenuation factor (i.e. rain over the radome)

    %minelev=1;

    %for ii=1:radar.ne  %Hago un loop sobre las elevaciones.
        %if( ii > 1 )
          %strhour=num2str(hour);
          %if(hour < 10);strhour=['0' strhour];end
          %strmin=num2str(min);
          %if(min < 10 );strmin=['0' strmin];end
          %%strsec=num2str(sec);
          %%if(sec < 10 );strsec=['0' strsec];end

          %fileout=[path '/ANGUIL2KM_LETKF_' date '_' strhour strmin '.dat'];
          %%I use the time corresponding to the first azimuth to decide if
          %%this scan is going to be included or not in the current group.
          %sec=sec + ( radar.time(1,ii) - radar.time(1,ii-1) );

        %end
        %if ( sec > 60 )

            %endelev= ii - 1;

            %tmp_radar=radar;
            %tmp_radar.elevation=radar.elevation(minelev:endelev);
            %tmp_radar.ne = ( endelev - minelev ) + 1;
            %tmp_radar.Z=radar.Z(:,:,minelev:endelev);
            %tmp_radar.longitude=radar.longitude(:,:,minelev:endelev);
            %tmp_radar.latitude =radar.latitude(:,:,minelev:endelev);
            %tmp_radar.Rh=radar.Rh(:,minelev:endelev);
            %tmp_radar.height=radar.height(:,minelev:endelev);


            %%fileout='tmp.bin';
            %fileout
            %[grid_ref , grid_count_ref , grid_dv , grid_count_dv , grid_az_ref , grid_el_ref , grid_ra_ref ]=radar_superobbing( tmp_radar , refc(:,:,minelev:endelev) , dvccc(:,:,minelev:endelev) , cart , fileout );

            %minelev=endelev + 1 ;

            %sec=sec-60;
            %min = min + 1;
            %if( min > 60 )
                %hour=hour + 1;
                %min=min-60;
            %end
            %if( hour > 23)
                %hour=0;
                %date=  datestr(datenum(date,'yyyymmdd')+1,'yyyymmdd');
            %end

        %end



    %end


    %Write the results in a binary data structure (For analysis and
    %forecast verification).
    %fileout=[FILEDIR '/ANGUIL2KM_LETKF_' date '_' strhour strmin '.rad'];
    %write_radar_data_seq(radar,refc,dvccc,attenuation_factor,qcflag,fileout,'b');


end
