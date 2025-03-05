#include "scalelib.h"
MODULE class_Ice_Shape
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! class_Ice_Shape
  ! version 4.0
  ! 08/06/2004
  !
  ! written by Tempei Hashino
  !
  ! All units are CGS.
  ! ergs = g cm^2/s^2
  !
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use scale_io
  use class_Thermo_Var
  use acc_amps
  use par_amps
  use mod_amps_const
  implicit none  
  private

  public :: make_Ice_Shape_Zero
  public :: ini_Ice_Shape_v4
  public :: get_total_sfc_area
  public :: get_perimeter
  public :: get_effect_area
  public :: get_circum_area
  public :: get_coef_ip
  public :: get_vip
  public :: get_vcs
  public :: cal_Ice_Shape_v4
  public :: cal_semiac_ip
  public :: cal_halfmaxdim_ip
  public :: diag_habit_v4

  public :: Ice_Shape

  TYPE Ice_Shape
     ! memory is aligned as declared below for use in common block
     sequence 
     ! --- variables necessary to discribe the 
     !     shape of ice hydrometeor ---
     !
     ! a + d = alen
     ! c
     ! 
     ! the a-axis length only on basal faces
     real(PS)    :: a
     ! the c-axis length
     real(PS)    :: c
     ! the dendrite arm
     real(PS)    :: d
     ! the rosette arm
     real(PS)    :: r
     ! the irregular crystal arm
     real(PS)    :: e


     ! the coordinate of center of gravity along a-axis 
     real(PS)    :: ag
     ! the coordinate of center of gravity along c-axis 
     real(PS)    :: cg
     ! the number of extra single ice crystals
     real(PS)    :: n_exice

     ! volume of a representative ice crystal
     real(PS)    :: V_ic

     ! volume of circumscribing spheroid of dry ice particle
     real(PS)    :: V_cs

     ! semi-major and minor axis lengths of dry ice particle
     real(PS)    :: semi_aip,semi_cip

     ! bulk sphere density of dry ice particle, and bulk crystal density of ice crystal
     real(PS)    :: den_ip,den_ic
     
     ! volume of circumscribing spheroid including water outside
     real(PS)    :: V_csw

!!c     ! mode for shape representation
!!c     ! 1 : spheroid
!!c     ! 2 : ice crystal model
!!c     integer     :: mode
     ! axis ratio of ice crystal, c/(a'+d)
     real(PS)                :: phi_ic
     ! axis ratio of ice crystal, d/a
     real(PS)                :: psi_ic
     ! axis ratio of ice crystal, r/a
     real(PS)                :: gam_ic
     ! axis ratio of ice crystal, e/a
     real(PS)                :: eta_ic

     ! axis ratio of circumscribing spheroid
     real(PS)                :: phi_cs     

     ! porocity or area ratio
     real(PS)                :: q_e


     ! habit of an ice crystal
     ! 1 : hexagonal plate
     ! 2 : broad-branch crystal
     ! 3 : columnar crystal
     ! 4 : rosette bullet
     ! 5 : irregular crystal (plate assemblage)
     integer     :: habit

     ! type of solid hydrometeor
     ! 1 : pristine crystal
     ! 2 : rimed crystal
     ! 3 : aggregate
     ! 4 : rimed aggregate
     ! 5 : graupel
     ! 6 : hail
     integer     :: sh_type

     ! ice particle model
     ! 1 : cylinder
     ! 2 : spheroid
     ! argument 1: ice crystal
     ! argument 2: ice particle
     integer,dimension(2)     :: is_mod
!tmp     integer,allocatable,dimension(:)     :: is_mod

     ! growth mode of ice crystals
     integer :: growth_mode

     ! indication of stochastic initialization of growth mode
     integer :: init_growth


  end type Ice_Shape

  ! number of lobes of the rosette
  real(PS), parameter         :: N_lob = 4.0

CONTAINS

  ! ************* Optional Constructor for all types **********************
!tmp  function make_Ice_Shape_zero () result (IS)
  subroutine make_Ice_Shape_zero(IS)
    type (Ice_Shape)  :: IS

    ! +++ calculate each component +++
    IS%a = 0.0_PS
    IS%d = 0.0_PS
    IS%c = 0.0_PS
    IS%r = 0.0_PS
    IS%e = 0.0_PS
    IS%ag = 0.0_PS
    IS%cg = 0.0_PS
    IS%n_exice = 0.0_PS

    ! +++ calculate the volume +++
    IS%V_ic = 0.0_PS
    ! circumscribing volume is calculated later.
    IS%V_cs = 0.0_PS

    IS%semi_aip = 0.0_PS
    IS%semi_cip = 0.0_PS
    IS%den_ip = 0.0_PS
    IS%den_ic = 0.0_PS
    IS%V_csw = 0.0_PS
    IS%q_e = 0.0_PS


    ! +++ categorize the ice hydrometeor +++
    IS%habit = 0
    IS%sh_type = 0

!tmp    call allocate_make_ice_shape
    IS%is_mod = 0


    IS%growth_mode=0
    IS%init_growth=0


!tmp  contains
!tmp    subroutine allocate_make_ice_shape
!tmp      integer :: item(1)
!tmp      item=0
!tmp      nullify(IS%is_mod)
!tmp      allocate( IS%is_mod(2), stat = item(1))
!tmp      if( any(item/=0) ) Stop 'Allocation failed at make_ice_shape'
!tmp    end subroutine allocate_make_ice_shape

!tmp  end function make_Ice_Shape_zero
  end subroutine make_Ice_Shape_zero

  subroutine undef_Ice_Shape(IS)
    type (Ice_Shape)  :: IS
    real(PS),parameter :: undef=999.9e+30
    integer,parameter :: iundef=999

    ! +++ calculate each component +++
    IS%a=undef
    IS%d=undef
    IS%c=undef
    IS%r=undef
    IS%e=undef
    IS%ag=undef
    IS%cg=undef
    IS%n_exice=undef
    IS%V_ic=undef
    IS%V_cs=undef
    IS%semi_aip=undef
    IS%semi_cip=undef
    IS%den_ip=undef
    IS%den_ic=undef
    IS%V_csw=undef
    IS%phi_ic=undef
    IS%psi_ic=undef
    IS%gam_ic=undef
    IS%eta_ic=undef
    IS%phi_cs=undef
    IS%q_e=undef
    IS%habit=iundef
    IS%sh_type=iundef
    IS%is_mod=iundef
    IS%growth_mode=iundef
    IS%init_growth=iundef
  end subroutine undef_Ice_Shape

!tmp  ! ***************** Public Deconstructor for a ice_shape type *****************
!tmp  subroutine delete_Ice_Shape (name)
!tmp    type (Ice_Shape),intent(inout) :: name
!tmp    integer :: item(1)
!tmp    item=0
!tmp    deallocate( name%is_mod, stat = item(1))
!tmp    nullify(name%is_mod)
!tmp    if( any(item/=0) ) Stop 'Deallocation failed at make_ice_shape'
!tmp    
!tmp  end subroutine delete_Ice_Shape

  subroutine cal_Ice_Shape ( level, th_var, IS, amass, alen, clen, em)
    use scale_prc, only: &
       PRC_abort
    use class_Thermo_Var
    ! level of complexity
    integer, intent(in)   :: level
    type (Thermo_Var), intent(in)  :: th_var
    type (Ice_Shape), intent(inout)  :: IS
    ! total mass of an ice crystal
    real(PS),pointer,dimension(:)  :: amass
    real(PS), dimension(3)    :: dum
    real(PS), intent(inout)             :: alen, clen
    real(PS) :: ratio,am4old,maxlen
    ! error message
    integer   :: em

    ! side plane crystal mass-length relation from Mitchel et al (1990)
    real(PS),parameter :: a_sdpl=0.00419,b_sdpl=2.3

    ! minimum (maximum) volume of circumscribing sphere corresponds to 1 um (10cm) radius
    real(PS),parameter :: V_csmin=4.18879020478639e-12,V_csmax=4.18879020478639e+3

!!c    write(*,'("1",5ES15.6)') IS%a,IS%c,IS%d,alen,clen

    ! +++ calculate the volume of ice crystal +++
    IS%V_ic = amass(imc)/den_i
!!c    write(*,'("mass(1:imc)",4ES15.6)') amass(imt),amass(imr),amass(ima),amass(imc)
!!c    if(max(alen,clen)>max(IS%e,IS%r)) then
    if( IS%V_ic < V_csmin ) then
!!c       write(*,*) " cal_ice_shape > volume is less than minimum:",IS%V_ic
!!c       write(*,'("mass(1:imc)",4ES15.6)') amass(imt),amass(imr),amass(ima),amass(imc)
       
       IS%V_ic=V_csmin
!!c       am4old=amass(imc)
       amass(imc)=IS%V_ic*den_i
       amass(imt)=max(amass(imt),amass(imc))
       if(amass(imr)+amass(ima)>1.0e-20_PS)then
          ratio=(amass(imt)-amass(imc))/(amass(imr)+amass(ima))
          amass(imr)=amass(imr)*ratio
          amass(ima)=amass(ima)*ratio
       else
          amass(imr)=0.0_PS
          amass(ima)=0.0_PS
       end if
!!c       if(level>=4) then
!!c          amass(imat)=amass(imat)*amass(imc)/am4old
!!c          amass(imas)=amass(imas)*amass(imc)/am4old
!!c          amass(imai)=amass(imai)*amass(imc)/am4old
!!c       end if
       alen=1.0e-4_PS
       clen=1.0e-4_PS
       IS%d=0.0_PS
       IS%r=0.0_PS
       IS%e=0.0_PS
       em=-1
    end if
!!c    if( amass(imt) == amass(imr)+amass(ima) ) then
!!c       write(*,'("1-1",5ES15.6)') amass(1:4)
!!c       ! --- this is the case of only rime or aggregate component ---
!!c       !     assume phi=1.0
!!c       amass(imc)=amass(imr)
!!c       amass(imr)=0.0_PS
!!c       IS%V_ic = amass(imc)/den_i
!!c       alen = ( IS%V_ic/(sq_three*3.0_PS))**(1.0/3.0)
!!c       clen = alen
!!c       IS%d = 0.0_PS
!!c    end if
    IS%a = alen-IS%d
    IS%c = clen

!!c    write(*,'("2",5ES15.6)') IS%a,IS%c,IS%d,alen,clen

    ! +++ validate the geometry of the ice crytal retrived from database +++
    call val_geo( level, IS, alen, clen, em)
    if( em >= 1 ) then
       return
    end if
!!c
!!c    alen = IS%a + IS%d
!!c    clen = IS%c
    IS%phi_ic = clen/alen
    IS%psi_ic = IS%d/alen
    IS%gam_ic = IS%r/alen
    IS%eta_ic = IS%e/alen

    ! +++ categorize the ice hydrometeor +++
    if( level == 1 ) then
       if( IS%phi_ic < 1.0_PS ) then
          IS%habit = 1
       else if( IS%phi_ic >= 1.0_PS ) then
          IS%habit = 3
       end if
    else if(level==2.or.level==4.or.level==6) then
       if( IS%phi_ic < 1.0_PS ) then
          if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
          else ; IS%habit = 2
          end if
       else
          IS%habit = 3       
       end if
    else if(level==3.or.level==5.or.level==7) then
       maxlen=max(alen,clen)
       if( IS%e>IS%r .and. IS%e>maxlen ) then
!!c       if( IS%e>IS%r .and. IS%e>0.0 ) then
          IS%habit = 5
!!c          write(*,*) "polycrystal!"
       elseif( IS%e<IS%r .and. IS%r>maxlen ) then
!!c       elseif( IS%e<IS%r .and. IS%r>0.0 ) then
          IS%habit = 4
!!c          write(*,*) "rosetta"
       else
          if( IS%phi_ic < 1.0_PS ) then
             if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
             else ; IS%habit = 2
             end if
          else
             IS%habit = 3       
          end if
       end if
    else
       LOG_ERROR("cal_Ice_shape",*) "This level for cal_Ice_Shape not defined:",level
       call PRC_abort
    end if
  end subroutine cal_Ice_Shape

  subroutine cal_Ice_Shape_v3 ( level, th_var, IS, amass, alen, clen, em)
    use scale_prc, only: &
       PRC_abort
    use class_Thermo_Var
    use mod_amps_utility, only: &
         get_len_s3, &
         get_len_s1, &
         get_len_c2a
    ! level of complexity
    integer, intent(in)   :: level
    type (Thermo_Var), intent(in)  :: th_var
    type (Ice_Shape), intent(inout)  :: IS
    ! total mass of an ice crystal
    real(PS),pointer,dimension(:)  :: amass
    real(PS), dimension(3)    :: dum
    real(PS), intent(inout)             :: alen, clen
    real(PS) :: ratio,am4old,maxlen
    ! error message
    integer   :: em

    ! minimum (maximum) volume of circumscribing sphere corresponds to 1 um (10cm) radius
    real(PS),parameter :: V_csmin=4.18879020478639e-12,V_csmax=4.18879020478639e+3
!!c    real(PS),parameter :: V_csmin=4.18879020478639e-12,V_csmax=4.18879020478639e+12

    ! factor
    real(PS),parameter :: fct1=0.1,fct2=0.9

    real(PS) :: vice1,vice2
    real(PS),parameter :: spx_p=0.25_PS ! sense test ros original spx_p=0.25, changed to 4

!!c    write(*,'("1",5ES15.6)') IS%a,IS%c,IS%d,alen,clen

    ! +++ calculate the volume of ice crystal +++
    vice1=amass(imc)/den_i
!!c    write(*,'("mass(1:4)",4ES15.6)') amass(imt),amass(imr),amass(ima),amass(imc)
!!c    if(max(alen,clen)>max(IS%e,IS%r)) then
    if(vice1<V_csmin) then
!!c       write(*,*) " cal_ice_shape > volume is less than minimum:",IS%V_ic
!!c       write(*,'("mass(1:4)",4ES15.6)') amass(imt),amass(imr),amass(ima),amass(imc)
       
       vice1=V_csmin
!!c       am4old=amass(imc)
       amass(imc)=vice1*den_i
       amass(imt)=max(amass(imt),amass(imc))
       if(amass(imr)+amass(ima)>1.0e-20_PS)then
          ratio=(amass(imt)-amass(imc))/(amass(imr)+amass(ima))
          amass(imr)=amass(imr)*ratio
          amass(ima)=amass(ima)*ratio
       else
          amass(imr)=0.0_PS
          amass(ima)=0.0_PS
       end if
       alen=1.0e-4_PS
       clen=1.0e-4_PS
       IS%d=0.0_PS
       IS%r=0.0_PS
       IS%e=0.0_PS
       IS%n_exice=0.0_PS
       IS%ag=0.0_PS
       IS%cg=0.0_PS
       em=-1
    end if
    IS%a = alen-IS%d
    IS%c = clen
    IS%phi_ic = clen/alen
    IS%psi_ic = IS%d/alen

!!c    write(*,'("2",5ES15.6)') IS%a,IS%c,IS%d,alen,clen

    ! +++ diagnose the habit and length +++
    if( level == 1 ) then
       if( IS%phi_ic < 1.0_PS ) then
          IS%habit = 1
       else if( IS%phi_ic >= 1.0_PS ) then
          IS%habit = 3
       end if
       vice2=coef4pi3*( alen**2.0+clen**2.0)**1.5
    else if(level==2.or.level==4.or.level==6) then
       if( IS%phi_ic < 1.0_PS ) then
          if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
          else ; IS%habit = 2
          end if
       else
          IS%habit = 3       
       end if
       vice2=coef4pi3*( alen**2.0+clen**2.0)**1.5
    else if(level==3.or.level==5.or.level==7) then
       if( IS%n_exice<0.5_PS) then
          ! single crystal
          if( IS%phi_ic < 1.0_PS ) then
             if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
             else ; IS%habit = 2
             end if
          else
             IS%habit = 3       
          end if
          IS%r=0.0_PS
          IS%e=0.0_PS

          vice2=coef4pi3*( alen**2.0+clen**2.0)**1.5
       else
          ! polycrystal
          if(IS%ag>=(IS%cg/clen+0.5_PS)*alen) then
             ! side planes or radiating assemblage of plates: S1 or P7a
             IS%habit = 5
             IS%e=get_len_s1(amass(imc))
             IS%r=0.0_PS

             vice2=coef4pi3*((1.0_ps+spx_p**2.0)**1.5)*IS%e**3.0
          elseif(IS%cg>=(IS%ag/alen+0.5_PS)*clen) then
             ! bullet rosette
             IS%habit = 4
             IS%r=get_len_c2a(amass(imc))
             IS%e=0.0_PS
             vice2=coef4pi3*((1.0_ps+spx_p**2.0)**1.5)*IS%r**3.0
          else
             ! combination of side planes, bullets, and columns: S3
             IS%habit = 6
             IS%e=get_len_s3(amass(imc))
             IS%r=IS%e
             vice2=coef4pi3*((1.0_ps+spx_p**2.0)**1.5)*IS%r**3.0
          end if

!!c          IS%r=sqrt(max(alen,IS%ag)**2.0+4.0_PS*max(clen,IS%cg)**2.0)
!!c          IS%e=sqrt(4.0_PS*max(alen,IS%ag)**2.0+max(clen,IS%cg)**2.0)
       end if

    else
       LOG_ERROR("cal_Ice_Shape_v3",*) "This level for cal_Ice_Shape not defined:",level
       call PRC_abort
    end if

    ! Obtain the volume of ice crystals
    IS%V_ic=max(vice1,vice2)
    IS%V_cs=min(V_csmax,max(IS%V_cs,IS%V_ic,amass(imt)/den_i))

    ! +++ validate the geometry of the ice crytal retrived from database +++
    call val_geo( level, IS, alen, clen, em)
    if( em >= 1 ) then
       return
    end if
!!c
!!c    alen = IS%a + IS%d
!!c    clen = IS%c
    IS%gam_ic = IS%r/alen
    IS%eta_ic = IS%e/alen

  end subroutine cal_Ice_Shape_v3

  subroutine cal_Ice_Shape_v4( level, th_var, IS, amass, alen, clen &
       ,vice1,em)
    use class_Thermo_Var
    use mod_amps_utility, only: get_len_s3,get_len_s1,get_len_c2a
    ! level of complexity
    integer, intent(in)   :: level
    type (Thermo_Var), intent(in)  :: th_var
    type (Ice_Shape), intent(inout)  :: IS
    ! total mass of an ice crystal
    real(PS),dimension(*)  :: amass
    real(PS), dimension(3)    :: dum
    real(PS), intent(inout)             :: alen, clen
    real(PS),intent(inout) :: vice1
    ! error message
    integer   :: em
    
    real(PS) :: ratio,am4old,maxlen

    ! minimum mass of hexagonal ice crystal with 1 um radius
    real(PS),parameter :: m_icmin=4.763209003e-12

    ! realistic bounds for length predictions
    real(ps),parameter :: min_hexlen=1.0e-4,max_hexlen=3.0

    ! factor
    real(PS),parameter :: fct1=0.1,fct2=0.9

    integer :: i_amass_gt0

!!c    write(*,'("1",5ES15.6)') IS%a,IS%c,IS%d,alen,clen

    ! +++ calculate the volume of ice crystal +++
    vice1=amass(imc)/den_i
!!c    write(*,'("mass(1:4)",4ES15.6)') amass(imt),amass(imr),amass(ima),amass(imc)
!!c    if(max(alen,clen)>max(IS%e,IS%r)) then
!!c    if(amass(imc)<m_icmin) then
    if(amass(imc)-m_icmin<-1.0e-4*m_icmin) then
!!c       write(*,*) " cal_ice_shape > volume is less than minimum:",IS%V_ic
!!c       write(*,'("mass(1:8)",8ES15.6)') amass(1:8)
       
       vice1=m_icmin/den_i
!!c       am4old=amass(imc)
       amass(imc)=vice1*den_i
       amass(imt)=max(amass(imt),amass(imc))

       i_amass_gt0=0.5_RP*(1.0_RP-sign(1.0_PS,1.0e-20_PS-(amass(imr)+amass(ima)+amass(imw))))
       ratio=(amass(imt)-amass(imc))/max(1.0e-20_RP,amass(imr)+amass(ima)+amass(imw))

       amass(imt)=real(i_amass_gt0,PS_KIND)*amass(imt) +&
                  (1.0_RP-real(i_amass_gt0,PS_KIND))*amass(imc)
       amass(imr)=real(i_amass_gt0,PS_KIND)*amass(imr)*ratio
       amass(ima)=real(i_amass_gt0,PS_KIND)*amass(ima)*ratio
       amass(imw)=real(i_amass_gt0,PS_KIND)*amass(imw)*ratio
       amass(imf)=0.0_PS

!org       if(amass(imr)+amass(ima)+amass(imw)>1.0e-20_PS)then
!org          ratio=(amass(imt)-amass(imc))/(amass(imr)+amass(ima)+amass(imw))
!org          amass(imr)=amass(imr)*ratio
!org          amass(ima)=amass(ima)*ratio
!org          amass(imw)=amass(imw)*ratio
!org          amass(imf)=0.0_PS
!org       else
!org          amass(imt)=amass(imc)
!org          amass(imr)=0.0_PS
!org          amass(ima)=0.0_PS
!org          amass(imw)=0.0_PS
!org          amass(imf)=0.0_PS
!org       end if

       alen=1.0e-4_PS
       clen=1.0e-4_PS
       IS%d=0.0_PS
       IS%r=0.0_PS
       IS%e=0.0_PS
       IS%n_exice=0.0_PS
       IS%ag=0.0_PS
       IS%cg=0.0_PS
       em=-1
    else
      alen=max(min_hexlen,min(max_hexlen,alen))
      clen=max(min_hexlen,min(max_hexlen,clen))
      IS%d=max(0.0_PS,min(0.9*alen,IS%d))
    end if
    IS%a = alen-IS%d
    IS%c = clen
!!!    IS%a = max(1.0e-4_PS,alen-IS%d)
!!!    IS%c = max(1.0e-4_PS,clen)
    IS%phi_ic = clen/alen
    IS%psi_ic = IS%d/alen

!!c    write(*,'("2",5ES15.6)') IS%a,IS%c,IS%d,alen,clen

    ! +++ diagnose the habit and length +++
!tmp    call diag_habit_v4(level,alen,clen,IS%d,IS%ag,IS%cg,IS%n_exice,&
!tmp         amass(imc),IS%r,IS%e,vice2,semi_a_i,semi_c_i,IS%habit)

!tmp    ! Obtain the volume of ice crystals
!tmp    IS%V_ic=max(vice1,vice2)
!tmp    IS%V_cs=min(V_csmax,max(IS%V_cs,IS%V_ic,amass(imt)/den_i))

!!c
!!c    alen = IS%a + IS%d
!!c    clen = IS%c
!tmp    IS%gam_ic = IS%r/alen
!tmp    IS%eta_ic = IS%e/alen

    ! +++ validate the geometry of the ice crytal retrived from database +++
!tmp    call val_geo( level, IS, alen, clen, em)
!tmp    if( em >= 1 ) then
!tmp       return
!tmp    end if
  end subroutine cal_Ice_Shape_v4

  subroutine ini_Ice_Shape_v4 (IS,alen,clen,semi_a_i,semi_c_i)
    type (Ice_Shape), intent(inout)  :: IS
    real(PS), intent(inout)             :: alen, clen,semi_a_i,semi_c_i

    IS%a=0.0_PS
    IS%c=0.0_PS
    IS%d=0.0_PS
    IS%r=0.0_PS
    IS%e=0.0_PS
    IS%ag=0.0_PS
    IS%cg=0.0_PS
    IS%n_exice=0.0_PS
    IS%V_ic=0.0_PS
    IS%V_cs=0.0_PS

    IS%phi_ic=0.0_PS
    IS%psi_ic=0.0_PS

    IS%habit=0


    alen=0.0_PS
    clen=0.0_PS
    semi_a_i=0.0_PS
    semi_c_i=0.0_PS
    IS%gam_ic=0.0_PS
    IS%eta_ic=0.0_PS

  end subroutine ini_Ice_Shape_v4


  subroutine Ini_Ice_Shape ( level, th_var, IS, tmass, alen, clen, em)
    use class_Thermo_Var
    ! level of complexity
    integer, intent(in)   :: level
    type (Thermo_Var), intent(in)  :: th_var
    type (Ice_Shape), intent(inout)  :: IS
    ! total mass of an ice crystal
    real(PS),intent(in)  :: tmass
    real(PS), dimension(3)    :: dum
    real(PS), intent(inout)             :: alen, clen
    ! error message
    integer   :: em

    if( tmass == 0.0 ) then
       ! +++ calculate each component +++
       alen = 0.0_PS
       clen = 0.0_PS
       IS%a = 0.0_PS
       IS%d = 0.0_PS
       IS%c = 0.0_PS
       ! +++ categorize the ice hydrometeor +++
       IS%habit = 0
       return
    end if

    dum(1)=(3.0**(-0.5))*(IS%V_ic)**(1.0/3.0)
    dum(2)=dum(1)
    dum(3)=0.0

    IS%a = dum(1)
    IS%c = dum(2)
    IS%d = dum(3)

    alen = IS%a + IS%d
    clen = IS%c
    IS%phi_ic = clen/alen
    IS%psi_ic = IS%d/alen

    ! +++ calculate the volume of ice crystal +++
!!c    IS%V_ic = tmass/den_i
    IS%V_ic = get_vol_pcd(alen,clen,IS%d)+get_vol_rss(IS%r)

    ! +++ categorize the ice hydrometeor +++
    if(level==1.or.level==2.or.level==4.or.level==6) then
       if( IS%phi_ic < 1.0_PS ) then
          IS%habit = 1
       else if( IS%phi_ic >= 1.0_PS ) then
          IS%habit = 3
       end if
    else if(level==3.or.level==5.or.level==7) then
       if( IS%phi_ic < 1.0_PS ) then
          if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
          else ; IS%habit = 2
          end if
       else
          IS%habit = 3       
       end if
    end if
  end subroutine Ini_Ice_Shape
    
  subroutine Ini_Ice_Shape_ros ( level, th_var, IS, amass, alen, clen, em)
    use class_Thermo_Var
    ! level of complexity
    integer, intent(in)   :: level
    type (Thermo_Var), intent(in)  :: th_var
    type (Ice_Shape), intent(inout)  :: IS
    ! total mass of an ice crystal
    real(PS),pointer,dimension(:)  :: amass
    real(PS), dimension(3)    :: dum
    real(PS), intent(inout)             :: alen, clen
    real(PS) :: maxlen
    ! error message
    integer   :: em



    IS%V_ic=5.196152423e-12_PS
    amass(imc)=IS%V_ic*den_i
    amass(imr)=0.0
    amass(ima)=0.0
    amass(imt)=amass(imc)


    ! +++ calculate each component +++
    alen = 5.0e-05_PS
    clen = 5.0e-05_PS


    IS%a = alen
    IS%d = 0.0_PS
    IS%c = clen


    IS%r=1.0e-04_PS
    IS%e=0.0_PS

    alen = IS%a + IS%d
    clen = IS%c
    IS%phi_ic = clen/alen
    IS%psi_ic = IS%d/alen


    IS%habit=4
    return


    ! +++ categorize the ice hydrometeor +++
    ! +++ categorize the ice hydrometeor +++
    if( level == 1 ) then
       if( IS%phi_ic < 1.0_PS ) then
          IS%habit = 1
       else if( IS%phi_ic >= 1.0_PS ) then
          IS%habit = 3
       end if
    else if( level == 2.or.level==4.or.level==6 ) then
       if( IS%phi_ic < 1.0_PS ) then
          if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
          else ; IS%habit = 2
          end if
       else
          IS%habit = 3       
       end if
    else if( level==3.or.level==5.or.level==7 ) then
       maxlen=max(alen,clen)
       if( IS%e>IS%r .and. IS%e>maxlen ) then
!!c       if( IS%e>IS%r .and. IS%e>0.0 ) then
          IS%habit = 5
!!c          write(*,*) "polycrystal!"
       elseif( IS%e<IS%r .and. IS%r>maxlen ) then
!!c       elseif( IS%e<IS%r .and. IS%r>0.0 ) then
          IS%habit = 4
!!c          write(*,*) "rosetta"
       else
          if( IS%phi_ic < 1.0_PS ) then
             if( IS%a > IS%d/2.0_PS ) then; IS%habit = 1
             else ; IS%habit = 2
             end if
          else
             IS%habit = 3       
          end if
       end if
    end if
  end subroutine Ini_Ice_Shape_ros
  subroutine Ini_Ice_Shape_hex ( level, th_var, IS, amass, alen, clen, em)
    use class_Thermo_Var
    ! level of complexity
    integer, intent(in)   :: level
    type (Thermo_Var), intent(in)  :: th_var
    type (Ice_Shape), intent(inout)  :: IS
    ! total mass of an ice crystal
    real(PS),pointer,dimension(:)  :: amass
    real(PS), dimension(3)    :: dum
    real(PS), intent(inout)             :: alen, clen
    real(PS) :: maxlen
    ! error message
    integer   :: em



    IS%V_ic=5.196152423e-12_PS
    amass(4)=IS%V_ic*den_i
    amass(2)=0.0
    amass(3)=0.0
    amass(1)=amass(4)


    ! +++ calculate each component +++
    alen = 1.0e-04_PS
    clen = 1.0e-04_PS


    IS%a = alen
    IS%d = 0.0_PS
    IS%c = clen


    IS%r=0.0_PS
    IS%e=0.0_PS

    alen = IS%a + IS%d
    clen = IS%c
    IS%phi_ic = clen/alen
    IS%psi_ic = IS%d/alen


    IS%habit=3
    return
  end subroutine Ini_Ice_Shape_hex

  function get_dendrite_area( IS, mode) result(A_d)
    ! +++ total area produced by dendritic arms which is projected from 
    !     the top of the basal plane +++
    type (Ice_Shape), intent(in)    :: IS
    integer                         :: mode
    real(PS)                           :: A_d
    integer :: i_mode_eq1

    i_mode_eq1=max(0.0,(-isign(1,mode-1)*max(0,mode-1)+1.0))*&
               max(0.0,(-isign(1,mode-1)*min(0,mode-1)+1.0))

    A_d= & ! case of mode == 1
        real(i_mode_eq1,PS_KIND)*0.0_PS +&
           ! case of mode == 2
        (1.0-real(i_mode_eq1,PS_KIND))*3.0*sq_three*IS%a*IS%d

!    if( mode == 1 ) then
!       A_d = 0.0_PS
!    else if( mode == 2) then
!       A_d = 3.0*sq_three*IS%a*IS%d
!    end if
  end function get_dendrite_area

  function get_basal_area( IS, mode) result(A_a)
    ! +++ area produced by the basal plane and projected from the 
    !     top of the basal plane +++
    type (Ice_Shape), intent(in)    :: IS
    integer                         :: mode
    real(PS)                           :: A_a
    integer :: i_mode_eq1

    i_mode_eq1=max(0.0,(-isign(1,mode-1)*max(0,mode-1)+1.0))*&
               max(0.0,(-isign(1,mode-1)*min(0,mode-1)+1.0))

    A_a= & ! case of mode == 1
        real(i_mode_eq1,PS_KIND)*PI*IS%a*IS%a +&
           ! case of mode == 2
        (1.0-real(i_mode_eq1,PS_KIND))*1.5*IS%a*IS%a*sq_three

!    if( mode == 1 ) then
!       A_a = PI*(IS%a**2.0)
!    else if( mode == 2) then
!       A_a = 1.5*(IS%a**2.0)*sq_three
!    end if
  end function get_basal_area

  function get_total_sfc_area( IS, alen, clen, mode) result(OMEGA)
    ! +++ gives total surface area +++
    type (Ice_Shape), intent(in)    :: IS
    integer                   :: mode
    real(PS)                  :: OMEGA
    real(PS)                  :: A_a, A_d
    real(PS), intent(in)      :: clen, alen
    real(PS)                  :: xx
    ! use Chiruta and Wang (2003)
!!c    real(PS),parameter        :: C_ros=1.5528*N_lob**0.7727
    real(PS),parameter        :: C_ros=4.532390052 ! sense test ros original 4.532390052, changed to 2.652903178
    integer :: i_mode_eq1,i_sh_le2,i_habit_eq4,i_ismod2_eq1,i_agtc,i_habit_le3 &
           ,i_aeqc,i_em1_le0
    real(PS) :: e1,e2,v_ros,v_cir,v_obl_1,v_obl_2,v_prl,v_hex

    i_mode_eq1=max(0.0_PS,(-isign(1,mode-1)*max(0,mode-1)+1.0_PS))*&
               max(0.0_PS,(-isign(1,mode-1)*min(0,mode-1)+1.0_PS))

    i_sh_le2=0.5_PS*(1.0_PS-real(isign(1,IS%sh_type-3),PS_KIND))

    i_habit_eq4=max(0.0_PS,(-isign(1,IS%habit-4)*max(0,IS%habit-4)+1.0_PS))*&
                max(0.0_PS,(-isign(1,IS%habit-4)*min(0,IS%habit-4)+1.0_PS))

    i_ismod2_eq1=max(0.0_PS,(-isign(1,IS%is_mod(2)-1)*max(0,IS%is_mod(2)-1)+1.0_PS))*&
                max(0.0_PS,(-isign(1,IS%is_mod(2)-1)*min(0,IS%is_mod(2)-1)+1.0_PS))

    i_agtc=0.5_PS*(1.0_PS-sign(1.0_PS,clen-alen))
    i_aeqc=max(0,floor(-sign(1.0_PS,alen-clen)*max(0.0_PS,alen-clen)+1.0_PS))*&
           max(0,floor(-sign(1.0_PS,alen-clen)*min(0.0_PS,alen-clen)+1.0_PS))

    i_habit_le3=0.5_PS*(1.0_PS-real(isign(1,IS%habit-4),PS_KIND))

    e1 = sqrt(1.0_PS - min(1.0_PS,(clen*clen)/(alen*alen)))
    e2 = sqrt(1.0_PS - min(1.0_PS,(alen*alen)/(clen*clen)))
    i_em1_le0=0.5_PS*(1.0_PS-sign(1.0_PS,e1-1.0_PS-1.0e-30_PS))

    A_a=1.5_PS*sq_three*(1.0_PS-IS%psi_ic)*alen*alen
    xx = 0.5_PS*alen*sqrt(1.0_PS+3.0_PS*IS%psi_ic*IS%psi_ic)

    v_ros=C_ros*IS%r*IS%r
    v_cir=4.0_PS*PI*alen*alen
    v_obl_1=2.0_PS*PI*alen*alen
    v_obl_2=2.0_PS*PI*alen*alen + PI*clen*clen* &
            log((1.0_PS+e1)/max(1.0e-30_PS,(1.0_PS-e1)))/max(1.0e-30_PS,e1)
    v_prl=2.0_PS*PI*alen*alen+2.0_PS*PI*alen*clen*asin(e2)/max(1.0e-30_PS,e2)
    v_hex=2.0_PS*A_a + 24.0_PS*clen*xx

    OMEGA= & ! case of mode == 1
        real(i_mode_eq1,PS_KIND)*&
        (real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND)*&
           v_ros+&
         (1.0_PS-real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND))*&
         ( real(i_aeqc,PS_KIND)*v_cir +&
           real(i_agtc,PS_KIND)* &
             (real(i_em1_le0,PS_KIND)*v_obl_1+(1.0_PS-real(i_em1_le0,PS_KIND))*v_obl_2) +&
           (1.0_PS-real(i_aeqc,PS_KIND))*(1.0_PS-real(i_agtc,PS_KIND))*   &
                   v_prl &
          )  &
        )+&
       ! case of mode==2
        (1.0_PS-real(i_mode_eq1,PS_KIND))*&
        ( &
          real(i_sh_le2,PS_KIND)*&
         (  real(i_habit_le3*i_ismod2_eq1,PS_KIND)*v_hex +&
            real(i_habit_eq4*i_ismod2_eq1,PS_KIND)*v_ros +&
        (1.0_PS-real(i_habit_le3*i_ismod2_eq1,PS_KIND))*(1.0_PS-real(i_habit_eq4*i_ismod2_eq1,PS_KIND))*&
         ( real(i_aeqc,PS_KIND)*v_cir +&
           real(i_agtc,PS_KIND)* &
                  (real(i_em1_le0,PS_KIND)*v_obl_1+(1.0_PS-real(i_em1_le0,PS_KIND))*v_obl_2) + &
           (1.0_PS-real(i_aeqc,PS_KIND))*(1.0_PS-real(i_agtc,PS_KIND))* &
                    v_prl &
          )) +&
        (1.0_PS-real(i_sh_le2,PS_KIND))*&
         ( real(i_aeqc,PS_KIND)*v_cir +&
           real(i_agtc,PS_KIND)* &
                  (real(i_em1_le0,PS_KIND)*v_obl_1+(1.0_PS-real(i_em1_le0,PS_KIND))*v_obl_2) + &
           (1.0_PS-real(i_aeqc,PS_KIND))*(1.0_PS-real(i_agtc,PS_KIND))* &
                    v_prl &
          ) &
        )

  end function get_total_sfc_area

  function get_perimeter( IS, alen, clen, mode) result(P)
    ! +++ perimeter normal to the flow +++
    type (Ice_Shape), intent(in)    :: IS
    real(PS), intent(in)               :: alen, clen
    integer                            :: mode
    real(PS)                  :: P
    real(PS)                  :: xx
    ! use Chiruta and Wang (2003)
    real(PS),parameter        :: C_ros=1.06655585 ! sens test ros original 2.1331117, changed to 1.06655585
    integer :: i_mode_eq1,i_sh_le2,i_habit_eq4,i_ismod2_eq1,i_agtc,i_habit_le3
    real(PS) :: v_ros,v_obl,v_prl,v_hex,v_col


    i_mode_eq1=max(0.0,(-isign(1,mode-1)*max(0,mode-1)+1.0))*&
               max(0.0,(-isign(1,mode-1)*min(0,mode-1)+1.0))

    i_sh_le2=0.5*(1.0-real(isign(1,IS%sh_type-3),PS_KIND))

    i_habit_eq4=max(0.0,(-isign(1,IS%habit-4)*max(0,IS%habit-4)+1.0))*&
                max(0.0,(-isign(1,IS%habit-4)*min(0,IS%habit-4)+1.0))

    i_ismod2_eq1=max(0.0,(-isign(1,IS%is_mod(2)-1)*max(0,IS%is_mod(2)-1)+1.0))*&
                max(0.0,(-isign(1,IS%is_mod(2)-1)*min(0,IS%is_mod(2)-1)+1.0))

    i_agtc=0.5*(1.0-sign(1.0_PS,clen-alen))

    i_habit_le3=0.5*(1.0-real(isign(1,IS%habit-4),PS_KIND))

    xx = sqrt( 0.25_PS*IS%a*IS%a + IS%d*IS%d - IS%a*IS%d/2.0_PS )

    v_ros=C_ros*IS%r
    v_obl=2.0_PS*PI*alen
    v_prl=PI*(3.0_PS*(alen+clen)-&
          sqrt((3.0_PS*alen + clen)*(alen + 3.0_PS*clen)))
    v_hex=12.0_PS*xx
    v_col=4.0_PS*(clen+alen)

    P = & ! case of mode == 1
      real(i_mode_eq1,PS_KIND)*&
        (  real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND)*v_ros+&
           (1.0-real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND))* &
            ( real(i_agtc,PS_KIND)*v_obl+(1.0-real(i_agtc,PS_KIND))*v_prl) ) &
        + &
      ! case of mode==2
      (1.0-real(i_mode_eq1,PS_KIND))*&
        ( real(i_sh_le2,PS_KIND)*&
           ( real(i_habit_le3*i_ismod2_eq1,PS_KIND)*&
              ( real(i_agtc,PS_KIND)*v_hex+(1.0-real(i_agtc,PS_KIND))*v_col) + &
             real(i_habit_eq4*i_ismod2_eq1,PS_KIND)*v_ros +&
             (1.0-real(i_habit_le3*i_ismod2_eq1,PS_KIND))*(1.0-real(i_habit_eq4*i_ismod2_eq1,PS_KIND))*&
              ( real(i_agtc,PS_KIND)*v_obl+(1.0-real(i_agtc,PS_KIND))*v_prl) ) +&
          (1.0-real(i_sh_le2,PS_KIND))* &
           ( real(i_agtc,PS_KIND)*v_obl+(1.0-real(i_agtc,PS_KIND))*v_prl) )

  end function get_perimeter
    
  function get_effect_area( IS,sh_type2,alen, clen, mode) result(A_e)
    ! +++ area projected from the top of the basal plane +++
    type (Ice_Shape), intent(in)    :: IS
    real(PS), intent(in)            :: alen, clen
    integer                         :: mode,sh_type2
    real(PS)                        :: A_e
    ! use Chiruta and Wang (2003) for 4 lobes
    real(PS),parameter        :: C_ros=0.899 ! sense test ros original C_ros=0.899, changed to 0.744
    integer :: i_mode_eq1,i_sh_le2,i_habit_eq4,i_ismod2_eq1,i_agtc,i_habit_le3
    real(PS) :: v_ros,v_obl,v_prl,v_hex,v_col

    i_mode_eq1=max(0.0,(-isign(1,mode-1)*max(0,mode-1)+1.0))*&
               max(0.0,(-isign(1,mode-1)*min(0,mode-1)+1.0))

    i_sh_le2=0.5*(1.0-real(isign(1,sh_type2-3),PS_KIND))

    i_habit_eq4=max(0.0,(-isign(1,IS%habit-4)*max(0,IS%habit-4)+1.0))*&
                max(0.0,(-isign(1,IS%habit-4)*min(0,IS%habit-4)+1.0))

    i_ismod2_eq1=max(0.0,(-isign(1,IS%is_mod(2)-1)*max(0,IS%is_mod(2)-1)+1.0))*&
                max(0.0,(-isign(1,IS%is_mod(2)-1)*min(0,IS%is_mod(2)-1)+1.0))

    i_agtc=0.5*(1.0-sign(1.0_PS,clen-alen))

    i_habit_le3=0.5*(1.0-real(isign(1,IS%habit-4),PS_KIND))

    v_ros=C_ros*IS%r*IS%r
    v_obl=PI*alen*alen
    v_prl=4.0_PS*clen*alen
    v_hex=1.5_PS*sq_three*(1.0_PS-IS%psi_ic)*alen*alen
    v_col=4.0_PS*clen*alen

    A_e = & ! case of mode == 1
        real(i_mode_eq1,PS_KIND)*&
          (&
          real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND)*v_ros +&
          (1.0_PS-real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND))*&
            (real(i_agtc,PS_KIND)*v_obl+(1.0_PS-real(i_agtc,PS_KIND))*v_prl) &
          )+&
       ! case of mode==2
        (1.0_PS-real(i_mode_eq1,PS_KIND))*&
          (&
          real(i_sh_le2,PS_KIND)*&
            (&
            real(i_habit_le3*i_ismod2_eq1,PS_KIND)*&
              (real(i_agtc,PS_KIND)*v_hex+&
              (1.0_PS-real(i_agtc,PS_KIND))*v_col)+&
            real(i_habit_eq4*i_ismod2_eq1,PS_KIND)*v_ros+&
            (1.0_PS-real(i_habit_le3*i_ismod2_eq1,PS_KIND))*(1.0_PS-real(i_habit_eq4*i_ismod2_eq1,PS_KIND))*&
              (real(i_agtc,PS_KIND)*v_obl+(1.0_PS-real(i_agtc,PS_KIND))*v_prl) &
            ) +&
          (1.0_PS-real(i_sh_le2,PS_KIND))*&
            (real(i_agtc,PS_KIND)*v_obl+(1.0_PS-real(i_agtc,PS_KIND))*v_prl)&
         )

  end function get_effect_area

  function get_circum_area( IS, sh_type2,alen, clen, mode) result(A_c)
    ! +++ area projected from the top of the basal plane +++
    type (Ice_Shape), intent(in)    :: IS
    real(PS), intent(in)               :: alen, clen
    integer                            :: mode,sh_type2
    real(PS)                           :: A_c
    integer :: i_mode_eq1,i_sh_le2,i_habit_eq4,i_ismod2_eq1,i_agtc
    real(PS) :: v_ros,v_obl,v_prl

    i_mode_eq1=max(0.0,(-isign(1,mode-1)*max(0,mode-1)+1.0))*&
               max(0.0,(-isign(1,mode-1)*min(0,mode-1)+1.0))
    i_sh_le2=0.5*(1.0-real(isign(1,sh_type2-3),PS_KIND))
    i_habit_eq4=max(0.0,(-isign(1,IS%habit-4)*max(0,IS%habit-4)+1.0))*&
                max(0.0,(-isign(1,IS%habit-4)*min(0,IS%habit-4)+1.0))
    i_ismod2_eq1=max(0.0,(-isign(1,IS%is_mod(2)-1)*max(0,IS%is_mod(2)-1)+1.0))*&
                max(0.0,(-isign(1,IS%is_mod(2)-1)*min(0,IS%is_mod(2)-1)+1.0))
    i_agtc=0.5*(1.0-sign(1.0_PS,clen-alen))

    v_ros=PI*IS%r*IS%r ! sense test ros original non multiplication, changed to /4
    v_obl=PI*alen*alen
    v_prl=PI*clen*alen

    A_c= &
        real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND)*v_ros +&
        (1.0_PS-real(i_sh_le2*i_habit_eq4*i_ismod2_eq1,PS_KIND))*&
          (real(i_agtc,PS_KIND)*v_obl+(1.0_PS-real(i_agtc,PS_KIND))*v_prl) 

  end function get_circum_area

!!c  function get_dm_d( IS, alen, clen, d_a, d_c, dm) result(dm_d)
!!c    ! +++ mass change in dendritic arms +++
!!c    type (Ice_Shape), intent(inout)    :: IS
!!c    real(PS), intent(in)               :: alen, clen, d_a, d_c
!!c    real(PS)                           :: dm_d, dm
!!c    ! assume that h_d = 2 c_len and delta d = delta a_len
!!c    dm_d = sq_three * IS%a * den_i * &
!!c         ( 2.0_PS * clen * d_a + IS%d * 2.0_PS * d_c )
!!c    if( dm_d > dm ) then
!!c       write(*,*) "Warning: dm_d > dm_t"
!!c       dm_d = dm
!!c       IS%h_d = ( dm_d/(sq_three*IS%a*den_i) - IS%d*2.0_PS*d_c)/d_a
!!c    end if
!!c  end function get_dm_d

  subroutine val_geo( level, IS, alen, clen, em)
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! validate the geometory of an ice crystal
    !   preserve the axis ratio from database
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)     :: level
    type (Ice_Shape), intent(inout)  :: IS
    real(PS), intent(inout) :: alen, clen
    ! error message
    integer   :: em
    real(PS)  :: dphi, dpsi

!!c    if( IS%a <= 1.0e-5_PS .and. IS%c <= 1.0e-5_PS ) then
    if( IS%a <= 0.0_PS .and. IS%c <= 0.0_PS ) then
!tmp       write(*,100) IS%a,IS%c,alen,clen
       em = 1
       return
100 format(" val_geo > the length for the mass components does not exists IS%a,IS%c,alen,clen:",4ES15.6)
    end if

    dphi = IS%c/alen
    dpsi = IS%d/alen

!!c    alen = ( IS%V_ic/(3.0_PS*sq_three*(1.0_PS-dpsi)*dphi ))**(1.0/3.0)

    IS%c = dphi*alen
    IS%d = dpsi*alen
    IS%a = alen-IS%d
    clen = IS%c
    if( alen == 0.0 .or. alen >= 1.0e+3) then
!tmp       write(*,'("vice,psi,phi,alen",4ES15.6)') IS%V_ic,dpsi,dphi,alen
       em=2
!!c       stop
    end if
  end subroutine val_geo

  subroutine cal_dadc_massc(is, dm_com, d_axis_len)
    type (Ice_Shape), intent(inout)  :: is
    ! change in a-axis and c-axis in the shifted bin
    real(PS), pointer, dimension(:)      :: d_axis_len
    ! mass change for a particle in the bin on
    ! 1: a axis
    ! 2: c axis
    ! 3: d axis
    real(PS), dimension(3)      :: dm_com

    d_axis_len(1) = dm_com(1)/(sq_three*den_i*is%c*(6.0_PS*is%a+2.0_PS*is%d))
    d_axis_len(2) = dm_com(2)/(sq_three*den_i*is%a*(3.0_PS*is%a+2.0_PS*is%d))
  end subroutine cal_dadc_massc

  function get_vol_pcd(a,c,d) result(out)
    real(PS),intent(in) :: a,c,d
    real(PS) :: psi,phi,out
    psi=d/a
    phi=c/a
    out=3.0_PS*sq_three*(a**3.0)*(1.0_PS-psi)*phi
  end function get_vol_pcd
  function get_vol_rss(r) result(out)
    real(PS),intent(in) :: r
    real(PS) :: out
    ! # of lobes of the rosette
    integer :: N
    N=4
    ! equation is from Chiruta and Wang (2003)
    out=0.3257*(real(N,PS_KIND)**0.5206)*(r**3.0)
  end function get_vol_rss

  subroutine diag_habit_v4(level,alen,clen,dlen,ag,cg,n_exice,&
       xlen_s1,xlen_c2a,xlen_s3,&
       rlen,elen,gam_ic,eta_ic,vice,semi_a_i,semi_c_i,habit)
    integer, intent(in) :: level
    real(PS), intent(in) :: alen,clen,dlen,ag,cg,n_exice
    real(PS), intent(in) :: xlen_s1,xlen_c2a,xlen_s3
    real(PS), intent(inout) :: rlen,elen,gam_ic,eta_ic,vice,semi_a_i,semi_c_i
    integer, intent(inout) :: habit

    real(PS),parameter :: spx_p=0.25_PS ! sense test ros original spx_p=0.25, changed to 4
    integer :: i_agtc,i_a2o3_gtd,i_ag_ge_cgclen,i_cg_ge_agalen


    i_agtc=0.5*(1.0-sign(1.0_PS,clen-alen))
    i_a2o3_gtd=0.5*(1.0-sign(1.0_PS,dlen-alen*2.0_PS/3.0_PS))
    i_ag_ge_cgclen=0.5*(1.0+sign(1.0_PS,ag-(cg/clen+0.5_PS)*alen))
    i_cg_ge_agalen=0.5*(1.0+sign(1.0_PS,cg-(ag/alen+0.5_PS)*clen))

!!!    i_ne_lt0p5=0.5*(1.0-sign(1.0,n_exice-0.5_PS))


    if( n_exice<0.5_PS) then
          ! single crystal
          habit=i_agtc*(i_a2o3_gtd*1+(1-i_a2o3_gtd)*2)+&
                (1-i_agtc)*3

          rlen=0.0_PS
          elen=0.0_PS

          vice=coef4pi3*( alen*alen+clen*clen)**1.5_PS

          semi_a_i=alen
          semi_c_i=clen

    else
          ! polycrystal
          habit=&
          ! side planes or radiating assemblage of plates: S1 or P7a
                i_ag_ge_cgclen*5 +&
          ! bullet rosette
                (1-i_ag_ge_cgclen)*i_cg_ge_agalen*4  +&
          ! combination of side planes, bullets, and columns: S3
                (1-i_ag_ge_cgclen)*(1-i_cg_ge_agalen)*6 

          elen=&
          ! side planes or radiating assemblage of plates: S1 or P7a
                real(i_ag_ge_cgclen,PS_KIND)*xlen_s1 +&
          ! bullet rosette
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*0.0_PS  +&
          ! combination of side planes, bullets, and columns: S3
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*(1.0_PS-real(i_cg_ge_agalen,PS_KIND))*&
                  xlen_s3

          rlen=&
          ! side planes or radiating assemblage of plates: S1 or P7a
                real(i_ag_ge_cgclen,PS_KIND)*0.0_PS +&
          ! bullet rosette
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*&
                  xlen_c2a  +&
          ! combination of side planes, bullets, and columns: S3
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*(1.0_PS-real(i_cg_ge_agalen,PS_KIND))*&
                  xlen_s3 

          semi_a_i=&
          ! side planes or radiating assemblage of plates: S1 or P7a
                real(i_ag_ge_cgclen,PS_KIND)*elen +&
          ! bullet rosette
                ! sense test ros
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*rlen +&
                !(1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*rlen/spx_p +&
          ! combination of side planes, bullets, and columns: S3
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*(1.0_PS-real(i_cg_ge_agalen,PS_KIND))*rlen

          semi_c_i=&
          ! side planes or radiating assemblage of plates: S1 or P7a
                real(i_ag_ge_cgclen,PS_KIND)*elen*spx_p +&
          ! bullet rosette
                ! sense test ros
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*rlen*spx_p +&
                !(1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*rlen +&
          ! combination of side planes, bullets, and columns: S3
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*(1.0_PS-real(i_cg_ge_agalen,PS_KIND))*rlen*spx_p

          vice=coef4pi3*((1.0_ps+spx_p**2.0_PS)**1.5_PS)

          vice=&
          ! side planes or radiating assemblage of plates: S1 or P7a
                real(i_ag_ge_cgclen,PS_KIND)*vice*elen*elen*elen +&
          ! bullet rosette
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*real(i_cg_ge_agalen,PS_KIND)*vice*rlen*rlen*rlen +&
          ! combination of side planes, bullets, and columns: S3
                (1.0_PS-real(i_ag_ge_cgclen,PS_KIND))*(1.0_PS-real(i_cg_ge_agalen,PS_KIND))*vice*rlen*rlen*rlen


!          if(ag>=(cg/clen+0.5_PS)*alen) then
!             ! side planes or radiating assemblage of plates: S1 or P7a
!             habit = 5
!             elen=get_len_s1(m_ice)
!             rlen=0.0_PS
!
!             vice=coef4pi3*((1.0_ps+spx_p**2.0)**1.5)*elen**3.0
!
!             semi_a_i=elen
!             semi_c_i=elen*spx_p
!
!          elseif(cg>=(ag/alen+0.5_PS)*clen) then
!             ! bullet rosette
!             habit = 4
!             rlen=get_len_c2a(m_ice)
!             elen=0.0_PS
!             vice=coef4pi3*((1.0_ps+spx_p**2.0)**1.5)*rlen**3.0
!
!             semi_a_i=rlen
!             semi_c_i=rlen*spx_p
!
!          else
!             ! combination of side planes, bullets, and columns: S3
!             habit = 6
!             elen=get_len_s3(m_ice)
!             rlen=elen
!             vice=coef4pi3*((1.0_ps+spx_p**2.0)**1.5)*rlen**3.0
!
!             semi_a_i=rlen
!             semi_c_i=rlen*spx_p
!
!          end if
!
    end if

    gam_ic = rlen/alen
    eta_ic = elen/alen


  end subroutine diag_habit_v4

  function get_vip(is_mod,phi_cs,semi_aip) result(v_ip)
    integer,intent(in) :: is_mod
    real(PS),intent(in) :: phi_cs,semi_aip
    real(PS) :: v_ip
    integer :: i_ismod_eq1

    i_ismod_eq1=max(0.0,(-isign(1,is_mod-1)*max(0,is_mod-1)+1.0))*&
                max(0.0,(-isign(1,is_mod-1)*min(0,is_mod-1)+1.0))

    v_ip= &
        ! is_mod==1
        !   assume the cylinder volume
        real(i_ismod_eq1,PS_KIND)*coef2p*phi_cs*semi_aip*semi_aip*semi_aip +&
        ! is_mod==2
        !   assume the spheroidal volume
        (1.0-real(i_ismod_eq1,PS_KIND))*coef4pi3*phi_cs*semi_aip*semi_aip*semi_aip
 

  end function get_vip

  function get_coef_ip(is_mod) result(coef)
    integer,intent(in) :: is_mod
    real(PS) :: coef
    if(is_mod==1) then
       ! assume the cylinder volume
       coef=coef2p
    else
       !if(is_mod==2) then
       ! assume the spheroidal volume
       coef=coef4pi3
    end if
  end function get_coef_ip
  function get_vcs(is_mod,phi_cs,semi_aip) result(v_cs)
    integer,intent(in) :: is_mod
    real(PS),intent(in) :: phi_cs,semi_aip
    real(PS) :: v_cs
    if(is_mod==1) then
       ! assume the cylinder volume
       v_cs=coef4pi3*semi_aip**3*(1.0_PS+phi_cs**2)**1.5
    else
       ! is_mod==2
       ! assume the spheroidal volume
       if(phi_cs<1.0_PS) then
          v_cs=coef4pi3*semi_aip**3
       else
          v_cs=coef4pi3*(semi_aip*phi_cs)**3
       end if
    end if
  end function get_vcs

!!c  function get_asr(is_mod,v_cs,semi_aip) result(phi_cs)
!!c    integer,intent(in) :: is_mod
!!c    real(PS),intent(in) :: v_cs,semi_aip
!!c    real(PS) :: phi_cs
!!c    if(is_mod==1) then
!!c       ! assume the cylinder volume
!!c       phi_cs=sqrt((v_cs/(coef4pi3*semi_aip**3.0))**(2.0/3.0)-1.0_PS)
!!c    elseif(is_mod==2) then
!!c       ! assume the spheroidal volume
!!c       if(i_cs<1.0_PS) then
!!c          v_cs=coef4pi3*semi_aip**3.0
!!c       else
!!c          v_cs=coef4pi3*(semi_aip*phi_cs)**3.0
!!c       end if
!!c    end if
!!c  end function get_vcs

  subroutine cal_semiac_ip(is_mod,phi_cs,v_cs,semi_aip,semi_cip)
    integer,intent(in) :: is_mod
    real(PS),intent(in) :: phi_cs,v_cs
    real(PS),intent(inout) :: semi_aip,semi_cip
    if(is_mod==1) then
       ! assume the cylinder volume
       semi_aip=(V_cs/coef4pi3)**(1.0/3.0)/sqrt(1.0_PS+phi_cs*phi_cs)
       semi_cip=semi_aip*phi_cs
    else   !! if(is_mod==2) then
       ! assume the spheroidal volume
       if(phi_cs<1.0) then
          semi_aip=(V_cs/coef4pi3)**(1.0/3.0)
          semi_cip=semi_aip*phi_cs
       else
          semi_cip=(V_cs/coef4pi3)**(1.0/3.0)
          semi_aip=semi_cip/phi_cs
       end if
    end if
  end subroutine cal_semiac_ip
  subroutine cal_halfmaxdim_ip(xlen,is_mod,semi_aip,semi_cip)
    real(PS),intent(inout) :: xlen
    integer,intent(in) :: is_mod
    real(PS),intent(in) :: semi_aip,semi_cip
    if(is_mod==1) then
      ! assume the cylinder volume
      xlen=sqrt(semi_aip*semi_aip+semi_cip*semi_cip)
    else !!!if(is_mod==2) then
      ! assume the spheroidal volume
      if(semi_aip>semi_cip) then
        xlen=semi_aip
      else
        xlen=semi_cip
      end if
    end if
  end subroutine cal_halfmaxdim_ip
  subroutine diag_is_mod(habit,sh_type,is_mod)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! determine the type of ice particle model
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer, intent(in)    ::  habit,sh_type
    integer,pointer,dimension(:) :: is_mod
    ! 1 : cylinder
    ! 2 : spheroid
    is_mod(1)=1
    if(sh_type==1) then
       is_mod(2)=1
    elseif(sh_type==2) then
       is_mod(2)=1
    elseif(sh_type==3) then
       is_mod(2)=1
    elseif(sh_type==4) then
       is_mod(2)=1
    elseif(sh_type==5) then
       is_mod(2)=2
    end if
  end subroutine diag_is_mod

END MODULE class_Ice_Shape
