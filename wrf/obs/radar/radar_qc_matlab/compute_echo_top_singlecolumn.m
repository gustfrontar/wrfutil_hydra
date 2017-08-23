

function [echo_top,echo_base,max_dbz,max_dbz_z,vertical_ref_grad]=compute_echo_top_singlecolumn(z,ref,mask,Nlevels,thresholds)

%Given a single reflectivity column compute:
%Echo_top : an Nlevels x Nthresholds array with the echo tops found in the
%column.
%Echo_base : as with echo top but for the echo bases.
%max_dbz : maximum dbz found in cloud regions defined by the thresholds and
%found at different levels.
%max_dbz_z : height of the maximum dbz for the different levels.

%Inputs: z, ref vectors with heights and reflectivities respectively.
%        mask : if 0 then the point is not used in the computations
%        (blocking or clutter detected by the mask, etc.).


%Start the computation

ref=squeeze(ref);
z=squeeze(z);
mask=squeeze(mask);

Nlevelstop=2; %How many levels with reflectivity below the threshold to declare an echo top.

echo_top=NaN(Nlevels,length(thresholds));
echo_base=NaN(Nlevels,length(thresholds));
max_dbz=NaN(Nlevels,length(thresholds));
max_dbz_z=NaN(Nlevels,length(thresholds));

vertical_ref_grad=NaN(length(z),1);

nan_profile=isnan( ref ) & (mask == 1);
ref(nan_profile)=-20;  %Asume a very low reflectivity value for these gates.

for it=1:length(thresholds)

base_detected=false;
top_detected=false;

base_count=0;

for ii=1:length(z)
 if(mask(ii) == 1 )   
   %Loop over the levels.
   
   %Look for an echo base
   
   if( ref(ii) > thresholds(it) && ~base_detected )
       %An echo base has been detected.
        
       if( base_count < Nlevels)
           
       base_detected=true;
       
       top_detected=false;
           
       base_count=base_count+1;
       
       echo_base(base_count,it)=z(ii);
       
       max_dbz(base_count,it)=ref(ii);
       
       max_dbz_z(base_count,it)=z(ii);

       end
       
   end
   
   %Look for an echo top.
  
   if ( ii > Nlevelstop )
      
   if ( sum ( ref(ii-Nlevelstop+1:ii) < thresholds(it) ) == Nlevelstop && ~top_detected && base_detected )
   %An echo top has been detected
       top_detected=true;
       
       base_detected=false;
       
       if( base_count <= Nlevels )
           
           echo_top(base_count,it)=z(ii-Nlevelstop);
           
       end
       
   end
   
   end
   
   %Echo top associated with top of the radar domain.
   
   if ( ii == length(z) && base_detected && ~top_detected )
    %Domain is over but echo top has not been found! :( 
    %Force echo top
      
       if( base_count <= Nlevels )
           echo_top( base_count ,it)=z(ii);
       end
       
   end
   
   %Compute max dbz
   if( base_detected && ~top_detected )
       %We are within a cloud or an echo region. Compute max dbz.
       if( ref(ii) > max_dbz(base_count,it) )
           
           max_dbz(base_count,it)=ref(ii);
           max_dbz_z(base_count,it)=z(ii);

       end
       
   end
   

   

 end %End over the if of the mask
    
end %End for loop over levels

end %End for loop over thresholds.

%Compute vertical profile of vertical ref gradient

for ii=1:length(z)-1;
    
    if(mask(ii) == 1)
       
        vertical_ref_grad(ii)=(ref(ii+1)-ref(ii))/(z(ii+1)-z(ii));
        
    end
    
    
end

vertical_ref_grad(end)=vertical_ref_grad(end-1);






