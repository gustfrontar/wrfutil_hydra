
function [attenuation_factor]=compute_atenuation_qc(radar,reflectivity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function estimates the attenuation percentaje due to metereological
%  echoes and also computes a correction.
%  Input:
%  radar an structure containing radar information.
%  reflectivity: a 3D array (possibly 4D) containing reflectivity data with
%  most of the ground clutter already removed (we will assume that the
%  reflectivity is associatedi with weather echoes only.
%  Output:
%  attenuation which is the attenuation factor A(r)  A. Berne and R. Uijlenhoet
%  2006.
%  This function is an approximate way to compute a quality control factor
%  for attenuation. The idea is to compute attenuation factor A(r) from the obtained
%  reflectivity. To avoid instability in the computaton of the attenuation
%  factor only the measured Z is used. This is an under estimation of the
%  attenuation but provides useful information for quality control purposes
%  (i.e. detection of strongly attenuated regions).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

% Define some parameters
%Perez and Zawasky 2003 a=0.000372, b=0.72 see Delrieu et al 2000 for other
%choices.

a=543;
b=1.36;
c=1.55e-3;
d=1.30;

alfa=(c^d)/(a^(d/b));
beta=(d/b);

beam_length=radar.range(2)-radar.range(1);

%Allocate variables
attenuation_factor=NaN(size(reflectivity));

Z=10.^(reflectivity/10);
Z(isnan(Z))=0;

%Check if we have information about the calibration error.
if( isfield(radar,'calibration_error') );
    if ( ~ isempty(radar.calibration_error) );
    calibration_error=radar.calibration_error;
    else
        calibration_error=1;
    end
else
    calibration_error=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative algorithm to compute attenuation. Based on the forward
% attenuation estimation of HB.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We iterate forward in the range direction.
attenuation_factor(:,1,:)=1;

%tmp=2*c*(Z/a).^(d/b);

for ir=1:radar.nr-1
    mean_k=squeeze(1e-3*(alfa*( (0.5)*(Z(:,ir,:)+Z(:,ir+1,:)) ).^beta ) );  %Compute mean k between ir and ir+1 (k is dbz/m);
    %attenuation_factor(:,ir+1,:)=sum(tmp(:,1:ir,:),2);
    attenuation_factor(:,ir+1,:)=squeeze(attenuation_factor(:,ir,:)).*exp(-0.46*mean_k*beam_length)*calibration_error;
end


time=toc;
display(['Atenuation was computed in ' num2str(time) ' seconds'])