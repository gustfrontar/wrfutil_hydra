!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Utility subroutines dealing super-droplet IDs
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-12 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2020-07-23 (S.Shima) [add] sdm_copy_selected_sd for sdm_dmpvar == 1?? and sdm_dmpvar == 2??
!! @li      2020-07-24 (S.Shima) [add] 'all', 'large', and 'selected' options to sdm_copy_selected_sd 
!! @li      2020-07-30 (S.Shima) [fix] sdm_copy_selected_sd 
!! @li      2020-10-30 (S.Shima) [fix] sdm_copy_selected_sd 
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module m_sdm_idutil
  use scale_precision
  use scale_io

  implicit none
  private
  public :: sdm_sort,sdm_getperm,sdm_copy_selected_sd

contains
  subroutine sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,    &
                         sort_tag0,fsort_tag,fsort_id,            &
                         sd_rand,sd_perm)

    ! Input variables
    integer, intent(in) :: freq_max ! maximum number of SD in each sd-grid
    integer, intent(in) :: ni_sdm   ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm   ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm   ! SDM model dimension in z direction
    integer, intent(in) :: sd_num   ! Number of super-droplets
    integer, intent(in) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)
    ! = sort_tag(n) - 1
    ! sort_tag(m) : accumulated number of
    !               super-droplets in each
    !               SDM-grid
    integer, intent(in) :: fsort_tag(0:freq_max+1)
    ! accumulated number with respect to
    ! the number contaiend super-droplets
    integer, intent(in) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
    ! super-droplets sorted by the number
    ! contaiend super-droplets
    ! Input and output variables
    real(DP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
    ! Output variables
    integer, intent(out) :: sd_perm(1:sd_num)    ! random permutations
    ! Work variables
    integer :: m, n, t, ss, tt         ! index
    !---------------------------------------------------------------

    ! Initialize for calculating random permutations
    
#ifdef QUICKDEBUG
    sd_perm(:) = -999
#endif
    !$omp parallel
    !$omp do private(t,m,tt)
    do t=fsort_tag(2),fsort_tag(freq_max+1)-1

       m = fsort_id(t)

       tt = sort_tag0(m) + 1
       sd_perm(tt) = 1

    end do
    !$omp end do

    ! Calculate random permutations

    do n=2,freq_max

       !$omp do private(t,m,tt,ss)
       do t=fsort_tag(n),fsort_tag(freq_max+1)-1

          m = fsort_id(t)

          tt = sort_tag0(m) + n
          !ORG        ss = sort_tag0(m) + ( int(sd_rand(tt)*n) + 1 )
          ss = sort_tag0(m) + ( int(sd_rand(tt)*n) + 1 )

          !### swap data ###!

          sd_perm(tt) = sd_perm(ss)
          sd_perm(ss) = n

       end do
       !$omp end do

    end do

    !$omp end parallel

    return
  end subroutine sdm_getperm
  !----------------------------------------------------------------------------
  subroutine sdm_sort(ni_sdm,nj_sdm,nk_sdm,                       &
                      sd_num,sd_n,sd_ri,sd_rj,sd_rk,                &
                      sort_id,sort_key,sort_freq,sort_tag,jdgtype)
    use scale_atmos_grid_cartesC_index, only: &
         IS,JS,KS
    use m_sdm_common, only: &
         VALID2INVALID,knum_sdm
    use gadg_algorithm, only: &
         gadg_count_sort
    ! Input variables
    integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
    integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
    integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    real(RP), intent(in) :: sd_ri(1:sd_num)  ! index[i/real] of super-droplets
    real(RP), intent(in) :: sd_rj(1:sd_num)  ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)  ! index[k/real] of super-droplets
    character(len=5), intent(in) :: jdgtype   ! flag for sorting
    ! Output variables
    integer, intent(out) :: sort_id(1:sd_num)   ! id that super-droplets sorted by sd-grids
    integer, intent(out) :: sort_key(1:sd_num)  ! sort key
    integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each sd-grid
    integer, intent(out) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)     ! accumulated number of super-droplets in each sd-grid
    integer :: adr                ! grid id
    integer :: max_key            ! total grid number
    integer :: i, j, k, n         ! index
    !--------------------------------------------------------------------

    ! Initialize

    max_key = ni_sdm * nj_sdm * knum_sdm + 1

    ! Sorting [step-1] -- initialize IKEY
    if( jdgtype .eq. 'valid' ) then

       !$omp parallel
       !$omp do private(n,i,j,k,adr) schedule(static,256) 
       do n=1,sd_num

          if( sd_rk(n)>VALID2INVALID ) then

             ! i \in {0,...,ni_sdm-1}, j \in {0,...,nj_sdm-1}, k \in {0,...,knum_sdm-1}
             i = floor(sd_ri(n)) - (IS-1) 
             j = floor(sd_rj(n)) - (JS-1)
             k = floor(sd_rk(n)) - (KS-1)

             adr = ni_sdm*nj_sdm*k + ni_sdm*j + i + 1

          else

             adr = max_key    !! invalid super-droplets

          end if

          sort_key(n) = adr

       end do
       !$omp end do
       !$omp end parallel

    else if( jdgtype .eq. 'multi' ) then

       !$omp parallel
       !$omp do private(n,i,j,k,adr) schedule(static,256) 
       do n=1,sd_num

          if( sd_rk(n)>VALID2INVALID .and. sd_n(n)>1 ) then

             i = floor(sd_ri(n)) - (IS-1) 
             j = floor(sd_rj(n)) - (JS-1)
             k = floor(sd_rk(n)) - (KS-1)

             adr = ni_sdm*nj_sdm*k + ni_sdm*j + i + 1

          else

             adr = max_key    !! invalid super-droplets

          end if

          sort_key(n) = adr

       end do
       !$omp end do
       !$omp end parallel

    end if

    
    ! Sorting [step-2] -- counting sort
    call gadg_count_sort( sort_key, 1, max_key,                       &
         sort_freq, sort_tag, sort_id )

    sort_tag(max_key+1) = sort_tag(max_key) + sort_freq(max_key)

    return
  end subroutine sdm_sort
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine sdm_copy_selected_sd(    &
       &                          KA, IA,IS,IE, JA,JS,JE, &
       &                          sd_num,    sd_numasl,    sd_n,    sd_x,    sd_y,    sd_ri,    sd_rj,    sd_rk,     &
       &                          sd_liqice,    sd_asl,    sd_r,    sdi,                                             &
       &                          sd_num_tmp,sd_numasl_tmp,sd_n_tmp,sd_x_tmp,sd_y_tmp,sd_ri_tmp,sd_rj_tmp,sd_rk_tmp, &
       &                          sd_liqice_tmp,sd_asl_tmp,sd_r_tmp,sdi_tmp,                                         &
       &                          TEMP0,ilist,sdtype)
    use scale_prc, only: &
         & PRC_abort
    use m_sdm_common, only: &
         & i2, sdicedef, sdm_cold, VALID2INVALID, STAT_LIQ, STAT_ICE, &
         & sdm_aslset, mass_amsul, ion_amsul, mass_nacl, ion_nacl, CurveF, ASL_FF 
    use m_sdm_coordtrans, only: &
         & sdm_x2ri, sdm_y2rj

    integer, intent(in) :: KA
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    integer,  intent(in)  :: sd_num      ! number of super-droplets
    integer,  intent(in)  :: sd_numasl   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    real(RP), intent(in)  :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in)  :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(in)  :: sd_rk(1:sd_num) ! face index-k(real) of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets

    integer,  intent(out) :: sd_num_tmp  ! number of super-droplets
    integer,  intent(out)  :: sd_numasl_tmp   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(DP), intent(out) :: sd_n_tmp(1:sd_num) ! multiplicity of super-droplets
    real(RP), intent(out)  :: sd_x_tmp(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(out)  :: sd_y_tmp(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri_tmp(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(out) :: sd_rj_tmp(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(out)  :: sd_rk_tmp(1:sd_num) ! face index-k(real) of super-droplets
    integer(i2), intent(out) :: sd_liqice_tmp(1:sd_num)
    real(RP), intent(out) :: sd_asl_tmp(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(out) :: sd_r_tmp(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(inout) :: sdi_tmp   ! ice phase super-droplets

    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! temperature [K]

    integer, intent(out) :: ilist(1:sd_num)   ! buffer for index list
    character(len=*), intent(in) :: sdtype

    real(RP) :: sd_thld_radi ! threshold radius [m]
    integer :: n,k,cnt,m,i,j,t,s
    real(RP):: sd_aslmw(1:22) ! Molecular mass of chemical material contained as water-soluble aerosol in super droplets (default+20)
    real(RP):: sd_aslion(1:22) ! Degree of ion dissociation of chemical material contained as water-soluble aerosol in super droplets (default+20)
    integer :: idx_nasl(1:22)  ! index for vactorization
    real(RP) :: dmask(1:22)  ! mask for vactorization
    real(RP):: coef_a, coef_b ! Coefficients of Kohler curve
    real(RP) :: t_sd      ! temperature of the grid contained the SD
    real(RP) :: ivt_sd    ! 1.d0 / t_sd
    real(RP) :: dtmp      ! temporary for interation  

    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)
    
    !### Copy the same aerosol chemical componets
    sd_numasl_tmp = sd_numasl

    !### Setup aerosol related parameters
    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then

          !### numasl=1 @ init : NaCl ###!

          sd_aslmw(1)  = mass_nacl
          sd_aslion(1) = ion_nacl

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslmw(1) = mass_amsul
          sd_aslmw(2) = mass_nacl
          sd_aslion(1) = ion_amsul
          sd_aslion(2) = ion_nacl

       end if

    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslmw(1) = mass_amsul
       sd_aslmw(2) = mass_nacl

       sd_aslion(1) = ion_amsul
       sd_aslion(2) = ion_nacl

       !! Must be a Bug. This cannot be simply commented out
       do n=1,20
       !            call getrname( id_sdm_aslmw  + (n-1), sd_aslmw(n+2)  )
       !            call getrname( id_sdm_aslion + (n-1), sd_aslion(n+2) )
       end do

    else if( abs(mod(sdm_aslset,10))==5 ) then

       !### numasl=1 @ init+rest : (NH4)HSO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    end if

    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do


    ! Get index list of the selected SDs
    cnt = 0
    if(sdtype == 'all') then
       do n=1,sd_num
          cnt = cnt + 1
          ilist(cnt) = n
       end do

    else if (sdtype == 'large') then
       sd_thld_radi = 1.0e-6_RP ! threshold radius [m]
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) cycle

          if( sd_liqice(n) == STAT_LIQ ) then
             if( sd_r(n)>sd_thld_radi ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if

          else if( sdm_cold .and. (sd_liqice(n) == STAT_ICE)) then
             if( (sdi%re(n)>sd_thld_radi) .or. (sdi%rp(n)>sd_thld_radi) ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if
          end if
       end do

    else if (sdtype == 'activated') then
       do n=1,sd_num
          if( sd_rk(n)<VALID2INVALID ) cycle

          if( sd_liqice(n) == STAT_LIQ ) then
             !! calculate the coefficient a of Kohler curve
             i = floor(sd_ri(n))+1
             j = floor(sd_rj(n))+1
             k = floor(sd_rk(n))+1

             t_sd  = TEMP0(k,i,j)

             ivt_sd = 1.0_RP / t_sd
             coef_a  = CurveF * ivt_sd

             !! calculate the coefficient b of Kohler curve
             coef_b = 0.0_RP

!OCL UNROLL('full'),NOSWP  
             do t=1,22

                s = idx_nasl(t)

                dtmp = sd_asl(n,s) * (real(sd_aslion(s),kind=RP)            &
                     / real(sd_aslmw(s),kind=RP))
                coef_b = coef_b + dmask(t) * dtmp
                
             end do

             coef_b = coef_b * ASL_FF

             !! calculate critical radius
             sd_thld_radi = sqrt(3.0_RP*coef_b/coef_a)

             if( sd_r(n)>sd_thld_radi ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if

          else if( sdm_cold .and. (sd_liqice(n) == STAT_ICE)) then
             sd_thld_radi = 1.0e-6_RP ! threshold radius [m]
             if( (sdi%re(n)>sd_thld_radi) .or. (sdi%rp(n)>sd_thld_radi) ) then
                cnt = cnt + 1
                ilist(cnt) = n
             end if
          end if
       end do

    else
       ! stop if unsupported sdtype option is specified
       LOG_ERROR("sdm_copy_selected_sd",*) "Unsupported sdtype option is specified"
       call PRC_abort

    end if
    sd_num_tmp = cnt

    ! Copy data of the selcted SDs
    if(sd_num_tmp /= 0) then
       !$omp parallel
       !$omp do simd private(n) schedule(static,256) 
       do m=1,sd_num_tmp
          n = ilist(m)
       
          sd_n_tmp(m)      = sd_n(n)
          sd_x_tmp(m)      = sd_x(n)
          sd_y_tmp(m)      = sd_y(n)
          sd_ri_tmp(m)     = sd_ri(n)
          sd_rj_tmp(m)     = sd_rj(n)
          sd_rk_tmp(m)     = sd_rk(n)
          sd_liqice_tmp(m) = sd_liqice(n)
          sd_r_tmp(m)      = sd_r(n)

       end do
       !$omp end do simd

       do k=1,sd_numasl_tmp
          !$omp do simd private(n) schedule(static,256) 
          do m=1,sd_num_tmp
             n = ilist(m)

             sd_asl_tmp(m,k) = sd_asl(n,k)

          end do
          !$omp end do simd
       end do

       if( sdm_cold ) then
          !$omp do simd private(n) schedule(static,256) 
          do m=1,sd_num_tmp
             n = ilist(m)

             sdi_tmp%re(m) = sdi%re(n)
             sdi_tmp%rp(m) = sdi%rp(n)
             sdi_tmp%rho(m) = sdi%rho(n)
             sdi_tmp%tf(m) = sdi%tf(n)
             sdi_tmp%mrime(m) = sdi%mrime(n)
             sdi_tmp%nmono(m) = sdi%nmono(n)

          end do
          !$omp end do simd
       end if
       !$omp end parallel
    end if

  end subroutine sdm_copy_selected_sd
end module m_sdm_idutil
