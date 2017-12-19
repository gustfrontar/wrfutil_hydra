import numpy as np
import datetime as dt
import os.path

default_undef_val=1.0e33

def read_data_direct(inputfilename,nx,ny,nz,dtypein,undef_in=default_undef_val,undef_out=default_undef_val,seq_acces=False):
#dtypein is the native input data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.

   if  os.path.exists(inputfilename) :

     f=open(inputfilename,'r')

     if ( seq_acces ) :

       field=np.ones((nx,ny,nz))*undef_out    

       for ii in range(0,nz) :
         nada=np.fromfile(f,dtype='>i4',count=1)
         tmpfield=np.fromfile(f,dtype=dtypein,count=nx*ny)
         nada=np.fromfile(f,dtype='>i4',count=1)
       
         field[:,:,ii]=np.reshape(tmpfield,(nx,ny))
       
     else :
       field=np.fromfile(f,dtype=dtypein,count=nx*ny*nz)

       field=np.reshape(field,(nz,nx,ny))  #Get a 3D-array
       field=field.transpose(1, 2, 0)      #Reorder array dimensions to (lon,lat,z)

     field[ abs(field) > undef_in ]=undef_out #Use the undef val of this module instead of the original undef value.

   else :

     print('Not found ',inputfilename)

     #If the file does not exist we will consider the entire data volume as missing data.
 
     field=np.zeros([nx,ny,nz])*undef_out


   return field

#def read_data_direct_woundef(inputfilename,nx,ny,nz,dtypein,undef_in=default_undef_val,undef_out=default_undef_val):
#dtypein is the native input data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.
#
#   if  os.path.exists(inputfilename) :
#   
#     f=open(inputfilename,'r')
#
#     field=np.fromfile(inputfilename,dtype=dtypein,count=nx*ny*nz)
#     field=np.reshape(field,(nz,nx,ny))  #Get a 3D-array
#     field=field.transpose(1, 2, 0)      #Reorder array dimensions to (lon,lat,z)
#
#   else :
#
#     field=np.zeros([nx,ny,nz])*undef_out
#
#   return field

def write_data_direct_woundef(inputfilename,field,dtypein):
#dtypein is the native input data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.
   tmp_shape=field.shape
   nx=tmp_shape[0]
   ny=tmp_shape[1]
   nz=tmp_shape[2]   

   #Reorder data to be consistent with input format.
   field=field.transpose(2,0,1)
   field=np.reshape(field,nx*ny*nz)

   #Write the data
   field.astype(dtypein).tofile(inputfilename)


def read_data_scale(basedir,expname,time,nx,ny,nz,undef_in=default_undef_val,undef_out=default_undef_val):

        datatypescale='>f4'

        data=dict()
        tmpvarname = 'U'
        filename=basedir  + tmpvarname + '_d01z-3d_' + time.strftime("%Y%m%d%H%M%S") + '.grd'

        data[tmpvarname ]=read_data_direct(filename,nx,ny,nz,datatypescale,undef_in=undef_in,undef_out=undef_out)

        tmpvarname = 'V'
        filename=basedir + tmpvarname + '_d01z-3d_' + time.strftime("%Y%m%d%H%M%S") + '.grd'

        data[tmpvarname ]=read_data_direct(filename,nx,ny,nz,datatypescale,undef_in=undef_in,undef_out=undef_out)

        tmpvarname = 'T'
        filename=basedir + tmpvarname + '_d01z-3d_' + time.strftime("%Y%m%d%H%M%S") + '.grd'

        data[tmpvarname ]=read_data_direct(filename,nx,ny,nz,datatypescale,undef_in=undef_in,undef_out=undef_out)

        tmpvarname = 'W'
        filename=basedir + tmpvarname + '_d01z-3d_' + time.strftime("%Y%m%d%H%M%S") + '.grd'

        data[tmpvarname ]=read_data_direct(filename,nx,ny,nz,datatypescale,undef_in=undef_in,undef_out=undef_out)

        tmpvarname = 'QV'
        filename=basedir + tmpvarname + '_d01z-3d_' + time.strftime("%Y%m%d%H%M%S") + '.grd'

        data[tmpvarname ]=read_data_direct(filename,nx,ny,nz,datatypescale,undef_in=undef_in,undef_out=undef_out)

        tmpvarname = 'QHYD'
        filename=basedir + tmpvarname + '_d01z-3d_' + time.strftime("%Y%m%d%H%M%S") + '.grd'

        data[tmpvarname ]=read_data_direct(filename,nx,ny,nz,datatypescale,undef_in=undef_in,undef_out=undef_out)

        return data


def read_data_scale_2(filename,nx,ny,nz,ctl_vars,ctl_inirecord,ctl_endrecord,dtypein='>f4',undef_in=default_undef_val,undef_out=default_undef_val):
#Default data type in is >f4 (double precission big endian)

#This function reads scale binary files generated with Guo-Yuan's scale python postprocessing tool.
#In these files all selected variables are included in the same binary file.
#nx,ny are the lon,lat dimensions nz is the total number of reccords vars*vert_levels 
#ctl_vars is a string list containig variable names.
#ctl_inirecord indicates the initial index for each variable (it size should be the same as ctl_vars)
#ctl_endrecord indicates the last index for each variable.
#A dictionary containing the different variables is returned. Each element of the dictionary is a masked array.

#Warning: Key order is not guaranteed so when a loop is performed over the keys, different key order may result
#in different excecutions of the code.

   my_data=dict()

   tmp_data=read_data_direct(filename,nx,ny,nz,dtypein,undef_in=undef_in,undef_out=undef_out)  #Read the data.

   #Loop over variables to create the dictionary. 
   ivar=0
   for my_var in ctl_vars   :

      my_data[my_var]=tmp_data[:,:,ctl_inirecord[ivar]:ctl_endrecord[ivar]+1]

      ivar=ivar+1

   return my_data






