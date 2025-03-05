Module maxdims
  use acc_amps
  implicit none
  public
  ! these parameters are set in init_AMPS subroutine
  !
  !
  ! max # of grid boxes in z, x, and y for the atmospheric part, 
  ! and z for the soil model.
  integer :: n1mx=125,n2mx=50,n3mx=2,n4mx=1
  integer :: nxpmax,nypmax,nzpmax,nzgmax

  integer :: mxln
  integer :: mxlh, mxlv
  !
  ! max # of grid systems
  integer,parameter :: maxgrds=2
  integer,parameter :: ntrgrds=maxgrds-1
  !
  ! # of total grib boxes that are passed to shared memory for 
  ! microphysical calculation. The best value depends on the machine.
  !integer,parameter :: lbnmx=200 ! for CASE_1
  !integer,parameter :: lensmx=lbnmx*330+50
  !
  ! # of soil types
  integer,parameter :: nstyp=12
  !
  ! 6D microphysical variables
  !    see the explanation at MASLWIKI.
  ! ***** This if for bulk micro par *****
  ! max # of parameters, bins, and categories for liquid 6d variable (qrp)
!  integer,parameter :: nprmx=4,nbrmx=1,ncrmx=2
  ! max # of parameters, bins, and categories for ice 6d variable (qip)
!  integer,parameter :: npimx=4,nbimx=1,ncimx=5
  ! max # of parameters, bins, and categories for aerosol 6d variable (qap)
!  integer,parameter :: npamx=4,nbamx=1,ncamx=2
  ! ***** This if for CLR *****
  ! max # of parameters, bins, and categories for liquid 6d variable (qrp)
!  integer,parameter :: nprmx=4,nbrmx=1,ncrmx=2
  ! max # of parameters, bins, and categories for ice 6d variable (qip)
!  integer,parameter :: npimx=4,nbimx=1,ncimx=3
  ! max # of parameters, bins, and categories for aerosol 6d variable (qap)
!  integer,parameter :: npamx=3,nbamx=1,ncamx=5
  ! ***** This if for AMPS *****
  ! max # of parameters, bins, and categories for liquid 6d variable (qrp)
  integer :: nprmx=6,nbrmx=40,ncrmx=1
  ! max # of parameters, bins, and categories for ice 6d variable (qip)
  integer :: npimx=18,nbimx=20,ncimx=1
  ! max # of parameters, bins, and categories for aerosol 6d variable (qap)
  integer :: npamx=5,nbamx=1,ncamx=4
  !
  ! max # of variable initialization files
  integer,parameter :: maxhfils=50
  !
  ! # of variable used in KUO cumulus parameterization
  !integer,parameter :: nkp=2*n1mx
  !
  ! critical (or minimum possible) distance between two grid interfaces in z direction.
  real(RP),parameter :: dzcrit=1.0_RP
  ! ---------------------------------------------------------------------
      
  ! ---------------------------------------------------------------------
  ! AMPS parameters for memory allocation.
  !   Other parameters are defined in par_micro.f90
  !
  ! number of AMPS objects 
  ! if running bulk microphysics parameterization, set LMAX=1.
  ! Otherwise, LMAX should be equal to lbin.
!  integer,parameter :: LMAX=1
  integer :: LMAX
  !
  ! max number of bins for each categories
  integer :: mxnbinr,mxnbini,mxnbina
  !
  ! max number of bins and bin boundaries for liq and ice hydrometeors
  integer :: mxnbin,mxnbinb
  !
  ! max number of mass components and volume components for ice particles
  integer,parameter :: mxnmasscomp=8,mxnvol=2
  !
  ! max number of microphysical tendencies, axis lengths, and non-mass component
  ! variables for ice particles.
  integer,parameter :: mxntend=12,mxnaxis=5,mxnnonmc=7
  !
  ! max number of mass components for liquid particles
  integer,parameter :: mxnmasscomp_r=3
  !
  ! ---------------------------------------------------------------------
end Module maxdims
