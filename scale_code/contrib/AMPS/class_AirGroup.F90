!OCL SERIAL
#include "scalelib.h"
MODULE class_AirGroup
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! class_AirGroup.f90
!
! version 1.0
!
! The object Air may consist of thermodynamics variables.
!
! written by Tempei Hashino
!
! All units are CGS.
! ergs = g cm^2/s^2
! 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  use scale_io
  use maxdims
  use class_Thermo_Var
  use acc_amps
  use mod_amps_const
  implicit none
  private

  public :: make_AirGroup_alloc
  public :: make_AirGroup_2
  public :: ini_AirGroup
  public :: update_AirGroup

  public :: AirGroup

  TYPE AirGroup
     ! memory is aligned as declared below for use in common block
     sequence 
!!c     private

     ! saturation vapor lookup table
     real(DS)  :: estbar(150),esitbar(111)
  
     ! thermo variables object has property of each gridbox
     ! first argument : grid number
     !type (Thermo_Var), dimension(LMAX)       :: TV
     type (Thermo_Var), allocatable       :: TV(:)

     ! number of grids that have microphyscis in it. 
     integer     :: L

     integer :: align_AirGroup ! CHIARUI

  endtype AirGroup


CONTAINS

  ! ************* Optional Constructor for a airgroup type **********************
  function make_AirGroup (dL,dest,desit,RV,DEN,PT,T,W) &
       result (a)
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer, intent(in) :: dL
    real(PS)    :: dest(150),desit(111)
    real(MP_KIND),dimension(dL),intent(in)  :: RV,DEN,PT,T,W
    type (AirGroup)  :: a

    ! optional status check
    integer  :: var_Status  
    integer  :: i, j

    allocate( a%TV(LMAX) )

    do i=1,dL
       if(DEN(i)<0.0)then
          LOG_ERROR("make_AirGroup",*) "class_air : Density is negative at",i
          do j=1,dL
             LOG_ERROR_CONT(*) j,DEN(j)
          end do
          call PRC_abort
       end if
    end do

    ! definition of total number of grids point with microphysics in it.
    a%L = dL

    ! initialization of lookup tables
    a%estbar=dest
    a%esitbar=desit

!tmp    ! +++ allocate memory to thermo variables +++
!tmp    allocate( a%TV(a%L), stat = var_Status)
!tmp    if( var_Status /= 0 ) stop "Memory not available for TV &
!tmp         in class_air"
    do i=1,a%L
       a%TV(i) = make_thermo_var3( RV(i),DEN(i),PT(i),T(i),W(i),a%estbar,a%esitbar )
    end do
  end function make_AirGroup

  function make_AirGroup_2 (dL,dest,desit,RV,DEN,PT,T,W,phase,RH) &
       result (a)
    implicit none
    integer, intent(in) :: dL,phase
    real(DS),intent(in)    :: dest(150),desit(111)
    real(PS),dimension(dL),intent(in)  :: RV,DEN,PT,T,W
    real(PS),intent(in) :: RH
    type (AirGroup)  :: a

    ! optional status check
    integer  :: var_Status  
    integer  :: i, j


    allocate( a%TV(LMAX) )

    ! definition of total number of grids point with microphysics in it.
    a%L = dL

    ! initialization of lookup tables
    a%estbar=dest
    a%esitbar=desit

!tmp    ! +++ allocate memory to thermo variables +++
!tmp    allocate( a%TV(a%L), stat = var_Status)
!tmp    if( var_Status /= 0 ) stop "Memory not available for TV &
!tmp         in class_air"
    do i=1,a%L
       a%TV(i) = make_thermo_var3_2(PT(i),T(i),W(i),a%estbar,a%esitbar,phase,RH )
    end do
  end function make_AirGroup_2
  function make_AirGroup_alloc (dL,dest,desit) &
       result (a)
    implicit none
    integer, intent(in) :: dL
    real(DS)    :: dest(150),desit(111)
    type (AirGroup)  :: a

    ! optional status check
    integer  :: var_Status  
    integer  :: i, j

    allocate( a%TV(LMAX) )

    ! definition of total number of grids point with microphysics in it.
    a%L = dL

    ! initialization of lookup tables
    a%estbar=dest
    a%esitbar=desit

    ! +++ allocate memory to thermo variables +++
!tmp    nullify(a%TV)
!tmp    allocate( a%TV(a%L), stat = var_Status)
!tmp    if( var_Status /= 0 ) stop "Memory not available for TV &
!tmp         in class_air"
  end function make_AirGroup_alloc

  subroutine ini_AirGroup(a,dL,RV,DEN,PT,T,W,rdsd,ihabit_gm_random)
    use scale_prc, only: &
       PRC_abort
    use mod_amps_utility, only: &
       cal_growth_mode_inl_vec, &
       random_genvar
    implicit none
    integer, intent(in) :: dL
    real(MP_KIND),dimension(DL),intent(in)  :: RV,DEN,PT,T,W
    type (AirGroup)  :: a
    type(random_genvar),intent(inout) :: rdsd
    ! random generaion: 1, max frequency: 0
    integer, intent(in)           :: ihabit_gm_random

    integer :: n
    integer,dimension(LMAX) :: ierror
    integer,dimension(LMAX) :: igm

    do n=1,dL
      ierror(n)=0
      if(DEN(n)<0.0)then
        ierror(n)=1
      endif
    enddo
    if(any(ierror(1:dL).gt.0)) then
       LOG_ERROR("ini_AirGroup",*) "class_air : Density is negative at",n
      do n=1,dL
         LOG_ERROR_CONT(*) n,DEN(n)
      end do
      call PRC_abort
    end if

    ! definition of total number of grids point with microphysics in it.
    a%L = dL

    do n=1,a%L
!       a%TV(n) = make_thermo_var3( RV(n),DEN(n),PT(n),T(n),W(n),a%estbar,a%esitbar )
      a%TV(n)%rv = rv(n)
      a%TV(n)%den = den(n)*1.0e-3_PS
      a%TV(n)%P = PT(n)*1.0e+1_PS
      a%TV(n)%T = T(n)
      a%TV(n)%W = W(n)*1.0e+2_PS

      a%TV(n)%D_v = get_diffusivity(a%TV(n)%P,a%TV(n)%T)
      a%TV(n)%k_a = get_thermal_conductivity(a%TV(n)%T)
      a%TV(n)%d_vis = get_dynamic_viscosity( a%TV(n)%T)
      a%TV(n)%sig_wa=get_sfc_tension(a%TV(n)%T)
      a%TV(n)%e=a%TV(n)%P*a%TV(n)%rv/(Rdvchiarui+a%TV(n)%rv)

      a%TV(n)%den_a=M_a*(a%TV(n)%P-a%TV(n)%e)/(R_u*a%TV(n)%T)

      a%TV(n)%dmassdt=0.0_PS

!!c      write(*,*) "ck rv in air:",i,a%TV(n)%rv

      a%TV(n)%e_sat(1) = get_sat_vapor_pres_lk(1, a%TV(n)%T, a%estbar, a%esitbar )
      a%TV(n)%e_sat(2) = get_sat_vapor_pres_lk(2, a%TV(n)%T, a%estbar, a%esitbar )

      a%TV(n)%s_v(1) = a%TV(n)%e/a%TV(n)%e_sat(1) - 1.0_PS
      a%TV(n)%GTP(1) = get_GTP( a%TV(n)%T, a%TV(n)%P, a%TV(n)%D_v, a%TV(n)%k_a, &
              a%TV(n)%e_sat(1), 1)
      a%TV(n)%rv_sat(1)=Rdvchiarui*a%TV(n)%e_sat(1)/max((a%TV(n)%P-a%TV(n)%e_sat(1)),a%TV(n)%e_sat(1))

      a%TV(n)%s_v(2) = a%TV(n)%e/a%TV(n)%e_sat(2) - 1.0_PS
      a%TV(n)%GTP(2) = get_GTP( a%TV(n)%T, a%TV(n)%P, a%TV(n)%D_v, a%TV(n)%k_a, &
              a%TV(n)%e_sat(2), 2)
      a%TV(n)%rv_sat(2)=Rdvchiarui*a%TV(n)%e_sat(2)/max((a%TV(n)%P-a%TV(n)%e_sat(2)),a%TV(n)%e_sat(2))

      if(a%TV(n)%rv_sat(1)<=a%TV(n)%rv) then
        a%TV(n)%s_v(1)=max(0.0_PS,a%TV(n)%s_v(1))
      else
        a%TV(n)%s_v(1)=min(0.0_PS,a%TV(n)%s_v(1))
      end if
      if(a%TV(n)%rv_sat(2)<=a%TV(n)%rv) then
        a%TV(n)%s_v(2)=max(0.0_PS,a%TV(n)%s_v(2))
      else
        a%TV(n)%s_v(2)=min(0.0_PS,a%TV(n)%s_v(2))
      end if

      a%TV(n)%dmassdt_v(1)=0.0_PS
      a%TV(n)%dmassdt_v(2)=0.0_PS
      a%TV(n)%dmassdt_v(3)=0.0_PS

      a%TV(n)%T_n=a%TV(n)%T
      a%TV(n)%T_m=a%TV(n)%T

      a%TV(n)%e_sat(2)=min(a%TV(n)%e_sat(1),a%TV(n)%e_sat(2))
      a%TV(n)%rv_sat(2)=min(a%TV(n)%rv_sat(1),a%TV(n)%rv_sat(2))

      a%TV(n)%s_v_n(1)=a%TV(n)%s_v(1)
      a%TV(n)%e_sat_n(1)=a%TV(n)%e_sat(1)
      a%TV(n)%s_v_n(2)=a%TV(n)%s_v(2)
      a%TV(n)%e_sat_n(2)=a%TV(n)%e_sat(2)
    enddo

!!c    write(*,*) "ck rdsd inair,gm",ihabit_gm_random
!!c    write(*,*) "ck rdsd inair0",rdsd
    call cal_growth_mode_inl_vec(igm,1,a%L,ihabit_gm_random &
                          ,a%TV(1:a%L)%T,a%TV(1:a%L)%s_v(2),rdsd)
!!c    write(*,*) "ck rdsd inair1",rdsd      

!!c    write(*,*) "rv1,rv2",a%TV(i)%rv,0.622*a%TV(i)%e/max((a%TV(i)%P-a%TV(i)%e),1.e-10)
    do n=1,a%L
      a%TV(n)%nuc_gmode=igm(n)
    enddo

  end subroutine ini_AirGroup

  subroutine update_AirGroup(a,RV,T,rdsd,ihabit_gm_random)
    use mod_amps_utility, only: cal_growth_mode_inl_vec,random_genvar
    implicit none
    real(PS),dimension(*),intent(in)  :: RV,T
    type (AirGroup)  :: a
    type(random_genvar),intent(inout) :: rdsd
    ! random generaion: 1, max frequency: 0
    integer, intent(in)           :: ihabit_gm_random

    integer :: n,in
    integer,dimension(LMAX) :: igm

    do n=1,a%L
      a%TV(n)%rv = rv(n)
      a%TV(n)%T = T(n)

      a%TV(n)%D_v = get_diffusivity(a%TV(n)%P,a%TV(n)%T)
      a%TV(n)%k_a = get_thermal_conductivity(a%TV(n)%T)
      a%TV(n)%d_vis = get_dynamic_viscosity( a%TV(n)%T)
      a%TV(n)%sig_wa=get_sfc_tension(a%TV(n)%T)
      a%TV(n)%e=a%TV(n)%P*a%TV(n)%rv/(Rdvchiarui+a%TV(n)%rv)

      a%TV(n)%den_a=M_a*(a%TV(n)%P-a%TV(n)%e)/(R_u*a%TV(n)%T)

      a%TV(n)%dmassdt=0.0_PS

      a%TV(n)%e_sat(1) = get_sat_vapor_pres_lk(1, a%TV(n)%T, a%estbar, a%esitbar )
      a%TV(n)%e_sat(2) = get_sat_vapor_pres_lk(2, a%TV(n)%T, a%estbar, a%esitbar )

      a%TV(n)%s_v(1) = min( a%TV(n)%e/a%TV(n)%e_sat(1) - 1.0_PS, 1000.0_PS )
      a%TV(n)%GTP(1) = get_GTP( a%TV(n)%T, a%TV(n)%P, a%TV(n)%D_v, a%TV(n)%k_a, &
              a%TV(n)%e_sat(1), 1)
      a%TV(n)%rv_sat(1)=Rdvchiarui*a%TV(n)%e_sat(1)/max((a%TV(n)%P-a%TV(n)%e_sat(1)),a%TV(n)%e_sat(1))

      a%TV(n)%s_v(2) = min( a%TV(n)%e/a%TV(n)%e_sat(2) - 1.0_PS, 1000.0_RP )
      a%TV(n)%GTP(2) = get_GTP( a%TV(n)%T, a%TV(n)%P, a%TV(n)%D_v, a%TV(n)%k_a, &
              a%TV(n)%e_sat(2), 2)
      a%TV(n)%rv_sat(2)=Rdvchiarui*a%TV(n)%e_sat(2)/max((a%TV(n)%P-a%TV(n)%e_sat(2)),a%TV(n)%e_sat(2))

      if(a%TV(n)%rv_sat(1)<=a%TV(n)%rv) then
        a%TV(n)%s_v(1)=max(0.0_PS,a%TV(n)%s_v(1))
      else
        a%TV(n)%s_v(1)=min(0.0_PS,a%TV(n)%s_v(1))
      end if
      if(a%TV(n)%rv_sat(2)<=a%TV(n)%rv) then
        a%TV(n)%s_v(2)=max(0.0_PS,a%TV(n)%s_v(2))
      else
        a%TV(n)%s_v(2)=min(0.0_PS,a%TV(n)%s_v(2))
      end if


      a%TV(n)%T_n=a%TV(n)%T
      a%TV(n)%T_m=a%TV(n)%T

      a%TV(n)%e_sat(2)=min(a%TV(n)%e_sat(1),a%TV(n)%e_sat(2))
      a%TV(n)%rv_sat(2)=min(a%TV(n)%rv_sat(1),a%TV(n)%rv_sat(2))

      a%TV(n)%dmassdt_v(1)=0.0_PS
      a%TV(n)%dmassdt_v(2)=0.0_PS
      a%TV(n)%dmassdt_v(3)=0.0_PS

      a%TV(n)%s_v_n(1)=a%TV(n)%s_v(1)
      a%TV(n)%e_sat_n(1)=a%TV(n)%e_sat(1)
      a%TV(n)%s_v_n(2)=a%TV(n)%s_v(2)
      a%TV(n)%e_sat_n(2)=a%TV(n)%e_sat(2)
    end do

    call cal_growth_mode_inl_vec(igm,1,a%L,ihabit_gm_random &
                          ,a%TV(1:a%L)%T,a%TV(1:a%L)%s_v(2),rdsd)

    do n=1,a%L
      a%TV(n)%nuc_gmode=igm(n)
    enddo

  end subroutine update_AirGroup

!tmp  subroutine delete_AirGroup (name)
!tmp    type (AirGroup), intent(inout) :: name
!tmp    ! check dealocate status
!tmp    integer                     :: ok
!tmp    integer                     :: i,j
!tmp
!tmp    deallocate (name%TV, stat = ok )
!tmp    if ( ok /= 0 ) stop "TV not deallocated in delete_AirGroup"
!tmp
!tmp  end subroutine delete_AirGroup

  subroutine renew_sat(ag,T,n,phase)
    type (AirGroup), intent(inout) :: ag
    ! temperature for which the saturation vapor is updated.
    real(PS),intent(in) :: T
    ! phase
    integer,intent(in) :: phase
    ! grid number
    integer,intent(in)  :: n
    integer  :: i
    select case (phase)
    case (1)
       i=1
       ag%TV(n)%e_sat(i) = get_sat_vapor_pres_lk(i, T, ag%estbar, ag%esitbar )
!!c       ag%TV(n)%rv_sat(i) = get_rv( T, ag%TV(n)%e_sat(i), ag%TV(n)%den)
!!c       RVIS(K)=.622*esi(k)/max((PT(K)-esi(k)),1.e-10)
       ag%TV(n)%rv_sat(i)=Rdvchiarui*ag%TV(n)%e_sat(i)/max((ag%TV(n)%P-ag%TV(n)%e_sat(i)),1.e-10_RP)
       ag%TV(n)%s_v(i) = ag%TV(n)%e/ag%TV(n)%e_sat(i) - 1.0_PS
       ag%TV(n)%GTP(i) = get_GTP( ag%TV(n)%T, ag%TV(n)%P, ag%TV(n)%D_v, ag%TV(n)%k_a, &
                         ag%TV(n)%e_sat(i), i)
    case (2)
       i=2
       ag%TV(n)%e_sat(i) = get_sat_vapor_pres_lk(i, T, ag%estbar, ag%esitbar )
       ag%TV(n)%rv_sat(i)=Rdvchiarui*ag%TV(n)%e_sat(i)/max((ag%TV(n)%P-ag%TV(n)%e_sat(i)),1.e-10_RP)
!!c       ag%TV(n)%rv_sat(i) = get_rv( T, ag%TV(n)%e_sat(i), ag%TV(n)%den)
       ag%TV(n)%s_v(i) = ag%TV(n)%e/ag%TV(n)%e_sat(i) - 1.0_PS
       ag%TV(n)%GTP(i) = get_GTP( ag%TV(n)%T, ag%TV(n)%P, ag%TV(n)%D_v, ag%TV(n)%k_a, &
                         ag%TV(n)%e_sat(i), i)
    case default
       do i = 1, 2
          ag%TV(n)%e_sat(i) = get_sat_vapor_pres_lk(i, T, ag%estbar, ag%esitbar )
          ag%TV(n)%rv_sat(i)=Rdvchiarui*ag%TV(n)%e_sat(i)/max((ag%TV(n)%P-ag%TV(n)%e_sat(i)),1.e-10_RP)
!!c          ag%TV(n)%rv_sat(i) = get_rv( T, ag%TV(n)%e_sat(i), ag%TV(n)%den)
          ag%TV(n)%s_v(i) = ag%TV(n)%e/ag%TV(n)%e_sat(i) - 1.0_PS
          ag%TV(n)%GTP(i) = get_GTP( ag%TV(n)%T, ag%TV(n)%P, ag%TV(n)%D_v, ag%TV(n)%k_a, &
                            ag%TV(n)%e_sat(i), i)
       end do
    end select

    ag%TV(n)%e_sat(2)=min(ag%TV(n)%e_sat(1),ag%TV(n)%e_sat(2))
    ag%TV(n)%rv_sat(2)=min(ag%TV(n)%rv_sat(1),ag%TV(n)%rv_sat(2))
  end subroutine renew_sat
  subroutine renew_sat_2(ag,T_old,esw,n,phase)
    type (AirGroup), intent(inout) :: ag
    ! temperature for which the saturation vapor is updated.
    real(PS),intent(in) :: T_old,esw
    ! phase
    integer,intent(in) :: phase
    ! grid number
    integer,intent(in)  :: n
    real(PS) :: T
    integer  :: i

    T= get_T_fesv_lk(1, T_old, esw, ag%estbar, ag%esitbar )
    select case (phase)
    case (1)
       i=1
       ag%TV(n)%e_sat(i) = get_sat_vapor_pres_lk(i, T, ag%estbar, ag%esitbar )
!!c       ag%TV(n)%rv_sat(i) = get_rv( T, ag%TV(n)%e_sat(i), ag%TV(n)%den)
!!c       RVIS(K)=.622*esi(k)/max((PT(K)-esi(k)),1.e-10)
       ag%TV(n)%rv_sat(i)=Rdvchiarui*ag%TV(n)%e_sat(i)/max((ag%TV(n)%P-ag%TV(n)%e_sat(i)),1.e-10_RP)
       ag%TV(n)%s_v(i) = ag%TV(n)%e/ag%TV(n)%e_sat(i) - 1.0_PS
       ag%TV(n)%GTP(i) = get_GTP( T, ag%TV(n)%P, ag%TV(n)%D_v, ag%TV(n)%k_a, &
                         ag%TV(n)%e_sat(i), i)
    case (2)
       i=2
       ag%TV(n)%e_sat(i) = get_sat_vapor_pres_lk(i, T, ag%estbar, ag%esitbar )
       ag%TV(n)%rv_sat(i)=Rdvchiarui*ag%TV(n)%e_sat(i)/max((ag%TV(n)%P-ag%TV(n)%e_sat(i)),1.e-10_RP)
!!c       ag%TV(n)%rv_sat(i) = get_rv( T, ag%TV(n)%e_sat(i), ag%TV(n)%den)
       ag%TV(n)%s_v(i) = ag%TV(n)%e/ag%TV(n)%e_sat(i) - 1.0_PS
       ag%TV(n)%GTP(i) = get_GTP( T, ag%TV(n)%P, ag%TV(n)%D_v, ag%TV(n)%k_a, &
                         ag%TV(n)%e_sat(i), i)
    case default
       do i = 1, 2
          ag%TV(n)%e_sat(i) = get_sat_vapor_pres_lk(i, T, ag%estbar, ag%esitbar )
          ag%TV(n)%rv_sat(i)=Rdvchiarui*ag%TV(n)%e_sat(i)/max((ag%TV(n)%P-ag%TV(n)%e_sat(i)),1.e-10_RP)
!!c          ag%TV(n)%rv_sat(i) = get_rv( T, ag%TV(n)%e_sat(i), ag%TV(n)%den)
          ag%TV(n)%s_v(i) = ag%TV(n)%e/ag%TV(n)%e_sat(i) - 1.0_PS
          ag%TV(n)%GTP(i) = get_GTP( T, ag%TV(n)%P, ag%TV(n)%D_v, ag%TV(n)%k_a, &
                            ag%TV(n)%e_sat(i), i)
       end do
    end select
!!c    if(ag%TV(n)%s_v(1)>0.0_PS) then
!!c       write(*,*) "T,T_old",T,T_old
!!c    end if
    ag%TV(n)%e_sat(2)=min(ag%TV(n)%e_sat(1),ag%TV(n)%e_sat(2))
    ag%TV(n)%rv_sat(2)=min(ag%TV(n)%rv_sat(1),ag%TV(n)%rv_sat(2))

  end subroutine renew_sat_2

  function get_satvp_lk(ag,phase,T) result(e_sat)
    type (AirGroup), intent(in) :: ag
    ! temperature for which the saturation vapor is updated.
    real(PS),intent(in) :: T
    ! grid number
    integer,intent(in)  :: phase
    real(PS) :: e_sat
    e_sat = get_sat_vapor_pres_lk(phase, T, ag%estbar, ag%esitbar )
  end function get_satvp_lk


END MODULE class_AirGroup

