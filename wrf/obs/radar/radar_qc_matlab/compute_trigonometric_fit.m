function [ alfa , beta , data_coverage , fit ] = compute_trigonometric_fit( f , tita  )
%COMPUTE_TRIGONOMETRIC_FIT ( radial velocity , azimuth , elevation )
% This function performs the trigonometric fit to the data over a single
% radar range to obtain an estimate of uo and vo that will be used in the
% construction of the vad profile.
% The algorithm is based on the derivative of vr with respect to tita as
% proposed by Gao and Drogenmeier 2003.

% Inputs:
% dvrdtita a vector with the derivatives of radial velocity with respect to
% the azimuth.
% tita: a vector with the corresponding azimuth values (in radians);
% elevation: the value of the current elevation angle.

% Note that dvrdtita could be an array of Naz , Nr, Nelev, in that case the
% fit is performed for each range and elevation, and the output values are
% matrices of size Nr x Nelev. Tita and elevation are still vectors with
% the different elevation angles.

% Output:
% The cos (alfa) and sin (beta) fit parameters
% Data_coverage : the percentage of points used to compute uo and vo at a
% certain elevation and range with respect to the total number of points at
% that elevation and range.

% The fit is performed using a trigonometric least squeares approach
% asuming a know frequency of 1;


%tic


if( size(f,1) == 1 );
    f=squeeze(f');
end
if( size(tita,1) == 1);
    tita=squeeze(tita');
end


[na , nr ,  ne ]= size(f);


alfa=NaN(nr,ne);
beta=NaN(nr,ne);
data_coverage=NaN(nr,ne);
fit=NaN(nr,ne);


for ir=1:nr
    for ie=1:ne
        tmp_f=f(:,ir,ie);
        tmp_index= ~isnan(tmp_f);

        if ( sum(tmp_index) > 5 ) 
        tmp_f=tmp_f( tmp_index );
        tmp_tita=tita( tmp_index );
        sintita=cos(tmp_tita);
        costita=sin(tmp_tita);

        M(1,1)=sum( costita.^2 );
        M(1,2)=sum( costita.*sintita );
        M(2,1)=M(1,2);
        M(2,2)=sum( sintita.^2 );
        
        B(1)=sum( costita.*tmp_f );
        B(2)=sum( sintita.*tmp_f );
        
        
        tmp_fit=(M^-1)*B';
        
        
        alfa(ir,ie)=tmp_fit(1);
        beta(ir,ie)=tmp_fit(2);
        
        %Compute the number of points used to compute uo and vo at each
        %range and elevation angle.
        data_coverage(ir,ie)=sum( tmp_index )/length(tita);
        
        tmp_data=alfa(ir,ie)*costita+beta(ir,ie)*sintita;
        fit(ir,ie)=var(tmp_data)/var(tmp_f); %Percentaje of explained variance.

        end
    end
end







%time=toc;
%display(['Trigonometric fitting was computed in ' num2str(time) ' seconds.'])




end

