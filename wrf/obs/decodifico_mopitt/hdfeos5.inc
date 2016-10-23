! * ==========================================================
! * File:       hdfeos5.inc
! * Author:     A. Muslimov, Emergent Technology Services, Inc.
! * Date:       March 11, 2002
! * Purpose:    Fortran header file for HDFEOS5 routines
! * ==========================================================
! * Contents: 
! *     HDFEOS5 data type tags
! *     HDFEOS5 logical constants
! *     HDFEOS5 file access tags
! *     HDFEOS5/GCTP library codes   
! *
! * ===========================================================
 
!          HDFEOS5  DATA TYPE TAGS 
!          =======================

           integer   HE5T_NATIVE_INT
           integer   HE5T_NATIVE_UINT
           integer   HE5T_NATIVE_SHORT
           integer   HE5T_NATIVE_USHORT
           integer   HE5T_NATIVE_SCHAR
           integer   HE5T_NATIVE_UCHAR
           integer   HE5T_NATIVE_LONG
           integer   HE5T_NATIVE_ULONG
           integer   HE5T_NATIVE_LLONG
           integer   HE5T_NATIVE_ULLONG
           integer   HE5T_NATIVE_FLOAT
           integer   HE5T_NATIVE_REAL
           integer   HE5T_NATIVE_DOUBLE
           integer   HE5T_NATIVE_LDOUBLE
           integer   HE5T_NATIVE_INT8
           integer   HE5T_NATIVE_UINT8
           integer   HE5T_NATIVE_INT16
           integer   HE5T_NATIVE_UINT16
           integer   HE5T_NATIVE_INT32
           integer   HE5T_NATIVE_UINT32
           integer   HE5T_NATIVE_INT64
           integer   HE5T_NATIVE_UINT64
           integer   HE5T_NATIVE_B8
           integer   HE5T_NATIVE_B16
           integer   HE5T_NATIVE_B32
           integer   HE5T_NATIVE_B64
           integer   HE5T_NATIVE_HSIZE
           integer   HE5T_NATIVE_HERR
           integer   HE5T_NATIVE_HBOOL

           integer   HE5T_STD_I8BE
           integer   HE5T_STD_I8LE
           integer   HE5T_STD_I16BE
           integer   HE5T_STD_I16LE
           integer   HE5T_STD_I32BE
           integer   HE5T_STD_I32LE
           integer   HE5T_STD_I64BE
           integer   HE5T_STD_I64LE
           integer   HE5T_STD_U8BE
           integer   HE5T_STD_U8LE
           integer   HE5T_STD_U16BE
           integer   HE5T_STD_U16LE
           integer   HE5T_STD_U32BE
           integer   HE5T_STD_U32LE
           integer   HE5T_STD_U64BE
           integer   HE5T_STD_U64LE
           integer   HE5T_STD_B8BE
           integer   HE5T_STD_B8LE
           integer   HE5T_STD_B16BE
           integer   HE5T_STD_B16LE
           integer   HE5T_STD_B32BE
           integer   HE5T_STD_B32LE
           integer   HE5T_STD_B64BE
           integer   HE5T_STD_B64LE

           integer   HE5T_IEEE_F32BE                      
           integer   HE5T_IEEE_F32LE                      
           integer   HE5T_IEEE_F64BE                      
           integer   HE5T_IEEE_F64LE                      
           integer   HE5T_NATIVE_CHAR
           integer   HE5T_CHARSTRING

           parameter(HE5T_NATIVE_INT      = 0)
           parameter(HE5T_NATIVE_UINT     = 1)
           parameter(HE5T_NATIVE_SHORT    = 2)
           parameter(HE5T_NATIVE_USHORT   = 3)
           parameter(HE5T_NATIVE_SCHAR    = 4)
           parameter(HE5T_NATIVE_UCHAR    = 5)
           parameter(HE5T_NATIVE_LONG     = 6)
           parameter(HE5T_NATIVE_ULONG    = 7)
           parameter(HE5T_NATIVE_LLONG    = 8)
           parameter(HE5T_NATIVE_ULLONG   = 9)
           parameter(HE5T_NATIVE_FLOAT    =10)
           parameter(HE5T_NATIVE_REAL     =10)
           parameter(HE5T_NATIVE_DOUBLE   =11)
           parameter(HE5T_NATIVE_LDOUBLE  =12)
           parameter(HE5T_NATIVE_INT8     =13)
           parameter(HE5T_NATIVE_UINT8    =14)
           parameter(HE5T_NATIVE_INT16    =15)
           parameter(HE5T_NATIVE_UINT16   =16)
           parameter(HE5T_NATIVE_INT32    =17)
           parameter(HE5T_NATIVE_UINT32   =18)
           parameter(HE5T_NATIVE_INT64    =19)
           parameter(HE5T_NATIVE_UINT64   =20)
           parameter(HE5T_NATIVE_B8       =21)
           parameter(HE5T_NATIVE_B16      =22)
           parameter(HE5T_NATIVE_B32      =23)
           parameter(HE5T_NATIVE_B64      =24)
           parameter(HE5T_NATIVE_HSIZE    =25)
           parameter(HE5T_NATIVE_HERR     =26)
           parameter(HE5T_NATIVE_HBOOL    =27)

           parameter(HE5T_STD_I8BE        =28)
           parameter(HE5T_STD_I8LE        =29)
           parameter(HE5T_STD_I16BE       =30)
           parameter(HE5T_STD_I16LE       =31)
           parameter(HE5T_STD_I32BE       =32)
           parameter(HE5T_STD_I32LE       =33)
           parameter(HE5T_STD_I64BE       =34)
           parameter(HE5T_STD_I64LE       =35)
           parameter(HE5T_STD_U8BE        =36)
           parameter(HE5T_STD_U8LE        =37)
           parameter(HE5T_STD_U16BE       =38)
           parameter(HE5T_STD_U16LE       =39)
           parameter(HE5T_STD_U32BE       =40)
           parameter(HE5T_STD_U32LE       =41)
           parameter(HE5T_STD_U64BE       =42)
           parameter(HE5T_STD_U64LE       =43)
           parameter(HE5T_STD_B8BE        =44)
           parameter(HE5T_STD_B8LE        =45)
           parameter(HE5T_STD_B16BE       =46)
           parameter(HE5T_STD_B16LE       =47)
           parameter(HE5T_STD_B32BE       =48)
           parameter(HE5T_STD_B32LE       =49)
           parameter(HE5T_STD_B64BE       =50)
           parameter(HE5T_STD_B64LE       =51)

           parameter(HE5T_IEEE_F32BE      =52)                
           parameter(HE5T_IEEE_F32LE      =53)                
           parameter(HE5T_IEEE_F64BE      =54)                
           parameter(HE5T_IEEE_F64LE      =55)                

           parameter(HE5T_NATIVE_CHAR     =56)                
           parameter(HE5T_CHARSTRING      =57)

           integer   HE5S_UNLIMITED_F
           parameter(HE5S_UNLIMITED_F = -1)


!          HDFEOS5 FILE ACCESS TAGS 
!          ========================
 
           integer   HE5F_ACC_RDWR
           integer   HE5F_ACC_RDONLY
           integer   HE5F_ACC_TRUNC

           parameter(HE5F_ACC_RDWR   = 100)
           parameter(HE5F_ACC_RDONLY = 101)
           parameter(HE5F_ACC_TRUNC  = 102)

!          TOOLKIT FILE ACCESS TAGS 
!          ========================
 
           integer   HDF4_ACC_RDONLY
           integer   HDF4_ACC_RDWR
           integer   HDF4_ACC_CREATE

           integer   HDF5_ACC_RDONLY
           integer   HDF5_ACC_RDWR
           integer   HDF5_ACC_CREATE

           parameter(HDF4_ACC_RDONLY = 1)
           parameter(HDF4_ACC_RDWR   = 3)
           parameter(HDF4_ACC_CREATE = 4)

           parameter(HDF5_ACC_RDONLY = 11)
           parameter(HDF5_ACC_RDWR   = 13)
           parameter(HDF5_ACC_CREATE = 14)


!          HDFEOS5  MERGE CODES 
!          ====================

           integer HE5_HDFE_NOMERGE
           integer HE5_HDFE_AUTOMERGE

           parameter(HE5_HDFE_NOMERGE   =0)
           parameter(HE5_HDFE_AUTOMERGE =1)
 
!          HDFEOS5 XXentries CODES 
!          =======================

           integer HE5_HDFE_NENTDIM
           integer HE5_HDFE_NENTMAP
           integer HE5_HDFE_NENTIMAP
           integer HE5_HDFE_NENTGFLD
           integer HE5_HDFE_NENTDFLD

           parameter(HE5_HDFE_NENTDIM  =0)
           parameter(HE5_HDFE_NENTMAP  =1)
           parameter(HE5_HDFE_NENTIMAP =2)
           parameter(HE5_HDFE_NENTGFLD =3)
           parameter(HE5_HDFE_NENTDFLD =4)

!          HDFEOS5 ANGLE CONVERSION CODES 
!          ==============================

           integer HE5_HDFE_RAD_DEG    
           integer HE5_HDFE_DEG_RAD    
           integer HE5_HDFE_DMS_DEG    
           integer HE5_HDFE_DEG_DMS    
           integer HE5_HDFE_RAD_DMS    
           integer HE5_HDFE_DMS_RAD    

           parameter(HE5_HDFE_RAD_DEG =0)
           parameter(HE5_HDFE_DEG_RAD =1)
           parameter(HE5_HDFE_DMS_DEG =2)
           parameter(HE5_HDFE_DEG_DMS =3)
           parameter(HE5_HDFE_RAD_DMS =4)
           parameter(HE5_HDFE_DMS_RAD =5)


!          HDFEOS5 SWATH SUBSET CODES 
!          ==========================
 
           integer HE5_HDFE_MIDPOINT 
           integer HE5_HDFE_ENDPOINT 
           integer HE5_HDFE_ANYPOINT 
           integer HE5_HDFE_INTERNAL 
           integer HE5_HDFE_EXTERNAL 
           integer HE5_HDFE_NOPREVSUB

           parameter(HE5_HDFE_MIDPOINT  =0) 
           parameter(HE5_HDFE_ENDPOINT  =1) 
           parameter(HE5_HDFE_ANYPOINT  =2) 
           parameter(HE5_HDFE_INTERNAL  =0) 
           parameter(HE5_HDFE_EXTERNAL  =1) 
           parameter(HE5_HDFE_NOPREVSUB =-1) 


!          HDFEOS5 GRID ORIGIN CODES 
!          =========================

           integer HE5_HDFE_GD_UL    
           integer HE5_HDFE_GD_UR    
           integer HE5_HDFE_GD_LL    
           integer HE5_HDFE_GD_LR    
 
           parameter(HE5_HDFE_GD_UL  =0) 
           parameter(HE5_HDFE_GD_UR  =1) 
           parameter(HE5_HDFE_GD_LL  =2) 
           parameter(HE5_HDFE_GD_LR  =3) 


!          HDFEOS5 PIXEL REGISTRATION CODES 
!          ================================

           integer HE5_HDFE_CENTER     
           integer HE5_HDFE_CORNER     
 
           parameter(HE5_HDFE_CENTER  =0)
           parameter(HE5_HDFE_CORNER  =1)

!          HDFEOS5 GCTP PROJECTION CODES 
!          =============================

           integer HE5_GCTP_GEO      
           integer HE5_GCTP_UTM      
           integer HE5_GCTP_SPCS     
           integer HE5_GCTP_ALBERS   
           integer HE5_GCTP_LAMCC     
           integer HE5_GCTP_MERCAT    
           integer HE5_GCTP_PS        
           integer HE5_GCTP_POLYC     
           integer HE5_GCTP_EQUIDC    
           integer HE5_GCTP_TM        
           integer HE5_GCTP_STEREO    
           integer HE5_GCTP_LAMAZ     
           integer HE5_GCTP_AZMEQD    
           integer HE5_GCTP_GNOMON    
           integer HE5_GCTP_ORTHO     
           integer HE5_GCTP_GVNSP     
           integer HE5_GCTP_SNSOID    
           integer HE5_GCTP_EQRECT    
           integer HE5_GCTP_MILLER    
           integer HE5_GCTP_VGRINT    
           integer HE5_GCTP_HOM       
           integer HE5_GCTP_ROBIN     
           integer HE5_GCTP_SOM       
           integer HE5_GCTP_ALASKA    
           integer HE5_GCTP_GOOD      
           integer HE5_GCTP_MOLL      
           integer HE5_GCTP_IMOLL     
           integer HE5_GCTP_HAMMER    
           integer HE5_GCTP_WAGIV     
           integer HE5_GCTP_WAGVII    
           integer HE5_GCTP_OBLEQA    
           integer HE5_GCTP_ISINUS    

           parameter(HE5_GCTP_GEO      =0)
           parameter(HE5_GCTP_UTM      =1)
           parameter(HE5_GCTP_SPCS     =2)
           parameter(HE5_GCTP_ALBERS   =3)
           parameter(HE5_GCTP_LAMCC    =4)
           parameter(HE5_GCTP_MERCAT   =5)
           parameter(HE5_GCTP_PS       =6)
           parameter(HE5_GCTP_POLYC    =7)
           parameter(HE5_GCTP_EQUIDC   =8)
           parameter(HE5_GCTP_TM       =9)
           parameter(HE5_GCTP_STEREO   =10)
           parameter(HE5_GCTP_LAMAZ    =11)
           parameter(HE5_GCTP_AZMEQD   =12)
           parameter(HE5_GCTP_GNOMON   =13)
           parameter(HE5_GCTP_ORTHO    =14)
           parameter(HE5_GCTP_GVNSP    =15)
           parameter(HE5_GCTP_SNSOID   =16)
           parameter(HE5_GCTP_EQRECT   =17)
           parameter(HE5_GCTP_MILLER   =18)
           parameter(HE5_GCTP_VGRINT   =19)
           parameter(HE5_GCTP_HOM      =20)
           parameter(HE5_GCTP_ROBIN    =21)
           parameter(HE5_GCTP_SOM      =22)
           parameter(HE5_GCTP_ALASKA   =23)
           parameter(HE5_GCTP_GOOD     =24)
           parameter(HE5_GCTP_MOLL     =25)
           parameter(HE5_GCTP_IMOLL    =26)
           parameter(HE5_GCTP_HAMMER   =27)
           parameter(HE5_GCTP_WAGIV    =28)
           parameter(HE5_GCTP_WAGVII   =29)
           parameter(HE5_GCTP_OBLEQA   =30)
           parameter(HE5_GCTP_ISINUS   =99)
 
 
!          HDFEOS5 TILIND/COMPRESSION CODES 
!          ================================

           integer HE5_HDFE_NOTILE       
           integer HE5_HDFE_TILE         
           integer HE5_HDFE_COMP_NONE    
           integer HE5_HDFE_COMP_RLE     
           integer HE5_HDFE_COMP_NBIT    
           integer HE5_HDFE_COMP_SKPHUFF 
           integer HE5_HDFE_COMP_DEFLATE 
           integer HE5_HDFE_COMP_SZIP_CHIP
           integer HE5_HDFE_COMP_SZIP_K13
           integer HE5_HDFE_COMP_SZIP_EC
           integer HE5_HDFE_COMP_SZIP_NN
           integer HE5_HDFE_COMP_SZIP_K13orEC
           integer HE5_HDFE_COMP_SZIP_K13orNN
           integer HE5_HDFE_COMP_SHUF_DEFLATE
           integer HE5_HDFE_COMP_SHUF_SZIP_CHIP
           integer HE5_HDFE_COMP_SHUF_SZIP_K13
           integer HE5_HDFE_COMP_SHUF_SZIP_EC
           integer HE5_HDFE_COMP_SHUF_SZIP_NN
           integer HE5_HDFE_COMP_SHUF_SZIP_K13orEC
           integer HE5_HDFE_COMP_SHUF_SZIP_K13orNN

           parameter(HE5_HDFE_NOTILE                 =0)
           parameter(HE5_HDFE_TILE                   =1)
           parameter(HE5_HDFE_COMP_NONE              =0)
           parameter(HE5_HDFE_COMP_RLE               =1)
           parameter(HE5_HDFE_COMP_NBIT              =2)
           parameter(HE5_HDFE_COMP_SKPHUFF           =3)
           parameter(HE5_HDFE_COMP_DEFLATE           =4)
           parameter(HE5_HDFE_COMP_SZIP_CHIP         =5)
           parameter(HE5_HDFE_COMP_SZIP_K13          =6)
           parameter(HE5_HDFE_COMP_SZIP_EC           =7)
           parameter(HE5_HDFE_COMP_SZIP_NN           =8)
           parameter(HE5_HDFE_COMP_SZIP_K13orEC      =9)
           parameter(HE5_HDFE_COMP_SZIP_K13orNN      =10)
           parameter(HE5_HDFE_COMP_SHUF_DEFLATE      =11)
           parameter(HE5_HDFE_COMP_SHUF_SZIP_CHIP    =12)
           parameter(HE5_HDFE_COMP_SHUF_SZIP_K13     =13)
           parameter(HE5_HDFE_COMP_SHUF_SZIP_EC      =14)
           parameter(HE5_HDFE_COMP_SHUF_SZIP_NN      =15)
           parameter(HE5_HDFE_COMP_SHUF_SZIP_K13orEC =16)
           parameter(HE5_HDFE_COMP_SHUF_SZIP_K13orNN =17)
 
!          HDFEOS5 GROUP CODES 
!          ===================

           integer HE5_HDFE_GEOGROUP        
           integer HE5_HDFE_DATAGROUP       
           integer HE5_HDFE_ATTRGROUP       
           integer HE5_HDFE_GRPATTRGROUP    
           integer HE5_HDFE_LOCATTRGROUP    
           integer HE5_HDFE_PROFGROUP       
           integer HE5_HDFE_PROFGRPATTRGROUP
           integer HE5_HDFE_GEOGRPATTRGROUP

           parameter(HE5_HDFE_GEOGROUP       =0)
           parameter(HE5_HDFE_DATAGROUP      =1)
           parameter(HE5_HDFE_ATTRGROUP      =2)
           parameter(HE5_HDFE_GRPATTRGROUP   =3)
           parameter(HE5_HDFE_LOCATTRGROUP   =4)
           parameter(HE5_HDFE_PROFGROUP      =5)
           parameter(HE5_HDFE_PROFGRPATTRGROUP =6)
           parameter(HE5_HDFE_GEOGRPATTRGROUP =7)
 















