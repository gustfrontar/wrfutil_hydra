import numpy as np
import numpy.ma as ma
import os

#Este modulo lee un ctl de grads para obtener parametros imporantes para la lectura de los datos

def read_ctl( filename , coding='utf-8' )  :

	fp=open(filename,'rb')
	ctl_ori=fp.read().decode(coding)
	ctl=dict()

	#Read some parameters and dimensions.

	ctl['template']=False
	ctl['big_endian']=False
	ctl['yrev']=False
	ctl['sequential']=False
	use_xdef=True
	for line in ctl_ori.split('\n'):
		if ( 'options' in line  or 'OPTIONS' in line ) and  ( 'template' in line or 'TEMPLATE' in line )  :
			ctl['template']=True
		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'big_endian' in line or 'BIG_ENDIAN' in line )  :
			ctl['big_endian']=True
		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'byteswapped' in line or 'BYTESWAPPED' in line )  :
			ctl['big_endian']=True

		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'sequential' in line or 'SEQUENTIAL' in line )  :
			ctl['sequential']=True
		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'yrev' in line or 'YREV' in line )  :
			ctl['yrev']=True

		if 'pdef' in line or 'PDEF' in line :
			ctl['nx']=float(line.split()[1])
			ctl['ny']=float(line.split()[2])
			ctl['dx']=float(line.split()[11])
			ctl['dy']=float(line.split()[12])
			use_xdef=False
		if ( 'xdef' in line or 'XDEF' in line ) and use_xdef :
			ctl['nx']=float(line.split()[1])  
		if ( 'ydef' in line or 'YDEF' in line ) and use_xdef :
			ctl['ny']=float(line.split()[1])  
		if 'zdef' in line or 'ZDEF' in line  :
			ctl['nz']=float(line.split()[1])
		if 'tdef' in line or 'TDEF' in line  :
			ctl['nt']=float(line.split()[1])
		if 'vars' in line or 'VARS' in line  and not 'ENDVARS' in line :
			ctl['nv']=float(line.split()[1])
		if 'undef' in line or 'UNDEF' in line :
			ctl['undef']=float(line.split()[1])

	#Read variable list
	ctl['var_list']=list()
	ctl['var_size']=list()
	ctl['var_desc']=list()

	var_section=False
	for line in ctl_ori.split('\n'):

		if ( 'ENDVARS' in line or 'endvars' in line ) :
			var_section=False
		if var_section                                :
			#Get the var name and store it in a var_list list	
			ctl['var_list'].append(line.split()[0])
			ctl['var_size'].append(line.split()[1])
			ctl['var_desc'].append(["".join(line.split()[2:])])
		if ( 'VARS' in line or 'vars' in line ) and not ( 'ENDVARS' in line or 'endvars' in line ) :
			var_section=True

	record=0
	ctl['ini_record']=list()
	ctl['end_record']=list()
	for var_size in ctl['var_size']   :
		ctl['ini_record'].append(record)
		if var_size == 0    :
			var_size = 1
		record=record + int(var_size) 
		ctl['end_record'].append(record-1)

	return ctl

def read_data_grads(filename,ctl,masked=False):

	my_data=dict()

	nx=int(ctl['nx'])
	ny=int(ctl['ny'])
	nz=int(np.max(ctl['end_record'])+1)

	undef=ctl['undef']

	if ctl['template']      :
		nt=1
	else                    :
		nt=int(ctl['nt'])

	if ctl['big_endian']    :
		dtypein='>f4'
	else                    :
		dtypein='f4'

	if ctl['sequential']    :
		sequential=True
	else                    :
		sequential=False

	tmp_data=read_data(filename,nx,ny,nz,nt,dtypein,undef,sequential)  #Read the data.

	#Loop over variables to create the dictionary. 
	for it in range(0,nt)      :
		ivar=0	
		for my_var in ctl['var_list']   :

			if it == 0        :
				nzvar=int(ctl['var_size'][ivar])
				if nzvar == 0    :
					nzvar=1
				my_data[my_var]=np.ones([ny,nx,nzvar,nt]).astype(np.float32)

			tmp_data_2=(tmp_data[:,:,ctl['ini_record'][ivar]:ctl['end_record'][ivar]+1,it])

			if masked                      :	
				rec=ctl['ini_record'][ivar]
				for iz in range(0,nzvar)  :
                        
					my_data[my_var][:,:,iz,it]=ma.masked_array( tmp_data[:,:,rec,it] , mask= tmp_data_2 == undef )
					rec=rec+1
 
			else                           :
				rec=ctl['ini_record'][ivar]
				for iz in range(0,nzvar)  :

					my_data[my_var][:,:,iz,it]=tmp_data[:,:,rec,it]
					rec=rec+1

			ivar=ivar+1

	tmp_data_2=None

	tmp_data=None
        


	return my_data



def read_data(inputfilename,nx,ny,nz,nt,dtypein,undef,seq_acces):
#dtypein is the native input data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.

	if  os.path.exists(inputfilename) :

		f=open(inputfilename,'r')

		field=np.ones((ny,nx,nz,nt)).astype(np.float32)*undef

		for it in range(0,nt) :
			for ii in range(0,nz) :
				if seq_acces :
					nada=np.fromfile(f,dtype='>i4',count=1)
				field[:,:,ii,it]= np.fromfile(f,dtype=dtypein,count=nx*ny).reshape(ny,nx) 
				if seq_acces :
					nada=np.fromfile(f,dtype='>i4',count=1)

	else :

		print('Not found ',inputfilename)

		#If the file does not exist we will consider the entire data volume as missing data.

		field=np.ones([ny,nx,nz,nt])*undef


	return field





	
        

