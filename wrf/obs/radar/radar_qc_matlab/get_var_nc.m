function out= get_var_nc(filename,varname)

%Obtenemos la cantidad de dimensiones, variables, atributos globales
ncid = netcdf.open(filename,'NC_NOWRITE');

[ndims nvars] = netcdf.inq(ncid);


%Obtengo el nombre de las dimensiones

for i=1:ndims
    
    
[dimname{i} dimlen(i)] = netcdf.inqDim(ncid,i-1);
    
end

%Obtengo el nombre de las variables

OFFSET=[];
SCALE=[];

for i=1:nvars
    
  [cvarname type tmpdims numatts] = netcdf.inqVar(ncid,i-1) ;
  
  if(strcmp(cvarname,varname)) 
  %Si esto se verifica, entonces esta es la variable que me interesa.
  numvar=i-1;
  
  tmpdata = double(netcdf.getVar(ncid,numvar));
  
      %Obtengo los nombres y valores de los atributos de la variable
      %seleccionada y los guardo en la estructura out.
      for j=1:numatts
          
        out.attributes(j).names = netcdf.inqAttName(ncid,numvar,j-1);
        out.attributes(j).values= netcdf.getAtt(ncid,numvar,out.attributes(j).names);

           %EN particular voy a usar los valores del add_offset y
           %scale_factor si estan presentes.
           
           if(strcmp(out.attributes(j).names,'add_offset'))
               OFFSET=double(out.attributes(j).values);
           end
           if(strcmp(out.attributes(j).names,'scale_factor'))
               SCALE=double(out.attributes(j).values);
           end
           
      end
      
      if(~(isempty(OFFSET) | isempty(SCALE)));
          
      display(['WARNING: SCALE FACTOR AND OFFSET ARE APPLIED TO DATA'])
       
      out.data=tmpdata*SCALE + OFFSET;
      
      else
      
      display(['WARNING: NO SCALE FACTOR AND OFFSET ATTRIBUTES HAS BEEN FOUND FOR THIS VARIBLE'])
      
      out.data=tmpdata;
  
      end
      
      
     %Ahora busco las dimensiones.
     
     for idim=1:length(tmpdims)
        out.dims(idim).name=dimname(tmpdims(idim)+1);
        out.dims(idim).len=dimlen(tmpdims(idim)+1);
        
        %Para cada dimension busco sus valores tambien.
        for n=1:nvars
          [tmpname] = netcdf.inqVar(ncid,n-1) ;
          if(strcmp(out.dims(idim).name,tmpname)) 
          %Si resulta que los valores que toma una dada dimension estan
          %presentes como una variable definida, entonces incluimos esa
          %informacion en la estructura.
            out.dims(idim).data = double(netcdf.getVar(ncid,n-1));
          end
        end
        
         
         
     end
      
      
  end    
end








