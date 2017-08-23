
function [attenuation_factor]=compute_atenuation(radar,reflectivity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function estimates the attenuation percentaje due to metereological
%  echoes and also computes a correction.
%  Input:
%  radar an structure containing radar information.
%  reflectivity: a 3D array (possibly 4D) containing reflectivity data with
%  most of the ground clutter already removed (we will assume that the
%  reflectivity is associatedi with weather echoes only.
%  Output:
%  attenuation which is the Path Integrated Attenuation (PIA) A. Berne and R. Uijlenhoet
%  2006.
%  correction is the correction (in dbz) that has to be added to each grid
%  point. 
%  To avoid the very well known instability of the forward attenuation
%  computation algorithm the algorithm is stopped when the attenuation
%  factor reaches a certain threshold. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

% Define some parameters
%Perez and Zawasky 2003 a=0.000372, b=0.72 see Delrieu et al 2000 for other
%choices.
a=0.000290; %A and b describes the relationship between reflectivity and attenuation factor.
b=0.72;
beam_length=radar.range(2)-radar.range(1);
instability_threshold=1/15;  %If attenuation factor is less than this then the algorithm stops.
%Delrieu et al 1999 suggests stop the correction when the reflectivity is
%incremented more than 10dbz for X band radars.

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

for ir=1:radar.nr-1
    %First: correct Z using the current value for the attenuation.
    Z(:,ir+1,:)=squeeze(Z(:,ir+1,:))./squeeze(attenuation_factor(:,ir,:));
    
    mean_k=squeeze(1e-3*(a/2)*(Z(:,ir,:).^b+Z(:,ir+1,:).^b));  %Compute mean k between ir and ir+1 (k is dbz/m);

    attenuation_factor(:,ir+1,:)=squeeze(attenuation_factor(:,ir,:)).*exp(-0.46*mean_k*beam_length)*calibration_error;
    
    %To avoid instability suppress attenuation computation if the beam
    %attenuation is larger than 0.5.
    %attenuation_factor(attenuation_factor < instability_threshold ) = NaN;
    
end

% Z(Z==0)=NaN;
% ref_cor=10*log10(Z);
% correction=ref_cor-reflectivity;

attenuation_factor=-10*log(attenuation_factor);  %Compute PIA


time=toc;
display(['Atenuation was computed in ' num2str(time) ' seconds'])