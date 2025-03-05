#include "scalelib.h"
!OCL SERIAL
module mod_amps_check
  use scale_io
  use maxdims
  use par_amps
  use com_amps, only: &
     fid_alog, &
     debug
  use mod_amps_const
  use class_AirGroup, only: &
     AirGroup
  use class_Group, only: &
     Group
  implicit none
  private

  public :: repair

contains

  subroutine repair(level,ag,gr,gs,ncat_a,ga,flagp_r,flagp_s,flagp_a,&
                    mes_rc, & !act_type,
                    irep_vap,irep_col,&
                    iupdate_gr,iupdate_gs,iupdate_ga,&
                    qtp,ID,JD,KD,from,itcall)
    use scale_prc, only: &
       PRC_abort
    !_______________________________________________________________________
    !     THIS ROUTINE CHECKS ALL COMPUTED CONVERSIONS TO SEE
    !     IF THEY WILL CAUSE A QUANTITY TO BECOME NEGATIVE OVER
    !     THE COURSE OF A TIMESTEP. IF SO, THE LOSSES ARE REDUCED.
    !     THIS IS AN ITERATIVE PROCESS BECAUSE THE LOSS FOR ONE
    !     CATEGORY IS A SOURCE FOR ANOTHER.
    !_______________________________________________________________________
    ! level of complexity
    integer, intent(in)           :: level
    !
    ! The modification constant is calculated for each group, not for
    ! each bin in a group.
    !
    type (AirGroup),intent(inout)    :: ag
    type (Group), intent(inout)         :: gr, gs

    integer,intent(in) :: ncat_a
    type (Group), dimension(ncat_a)         :: ga
    ! message from reality-check
!tmp    integer,pointer,dimension(:)  :: mes_rc
    integer,dimension(*)   :: mes_rc
    !integer, intent(in)           :: act_type
    real(MP_KIND) :: qtp(*)
    integer :: ID(*),JD(*),KD(*)
    character*(*) :: from
    integer, intent(in)           :: itcall
    ! repair flag for vapor and collision processes
    logical,intent(in) :: irep_vap,irep_col

    integer,dimension(*),intent(in) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,*),intent(in) :: iupdate_ga
    ! flag for prediction
    integer, intent(in)                   :: flagp_a, flagp_r, flagp_s

    !
    ! local space
    !
    ! accumulated (multiplied over the iteration) modification constant
    ! for rain group
    real (PS), dimension(mxnbin,mxntend,LMAX)    :: acc_mod_r
    ! for solid hydrometeor group
    real (PS), dimension(mxnbin,mxntend,LMAX)    :: acc_mod_s
    ! for aerosol group
    real (PS), dimension(mxnbina,mxntend,ncamx,LMAX)    :: acc_mod_a
    ! mark of modification
    integer,dimension(LMAX)   :: mark_mod
    ! mark for calculation of other tendencies
    integer,dimension(LMAX)   :: mark_all
    integer,dimension(LMAX)   :: mark_mass_r, mark_mass_s
    integer,dimension(LMAX)   :: mark_con_r, mark_con_s
    integer,dimension(LMAX)   :: mark_v_r, mark_v_s
    integer,dimension(ncamx,LMAX)    :: mark_mass_a,mark_con_a,mark_v_a
    ! maximum number of iteration
    integer                   :: NITER, NITER_TOTAL

    real(PS),dimension(LMAX)    :: rr_bef,rs_bef,cr_bef,rr,rs,cr,NDXDT,NDXDT2
    real(PS),dimension(LMAX)    :: cs_bef,cs
    real(PS) :: nxs,nxr
!!!    real(PS),dimension(mxnbini,mxntend) :: temp_dmdt_s,temp_dmdtap_s
!!!    real(PS),dimension(mxnbinr,mxntend) :: temp_dmdt_r,temp_dcdt_r

    real(PS),dimension(LMAX) :: udepr,uactr,udepi,uacti,udepr2,uactr2,udepi2,uacti2,&
          garm,lsmlt,garm2,lsmlt2,garml,garml2,lsmltl,lsmltl2,&
          udepar,udepai,udepar2,udepai2, &
          udhfi,udhfi2,&
          uactar,uactar2,uimmar,uimmar2,uimmai,uimmai2
    real(PS),dimension(ncamx,LMAX) :: udepa,udepa2
    real(PS),dimension(ncamx,LMAX) :: uacta,uacta2

!    integer,dimension(LMAX*mxnbin) :: ierror1
    integer,dimension(LMAX) :: icond1
    !
    !integer                   :: var_Status
    integer                   :: i,m,n, nvol,j,k,ica,itr1,itr2!, ngrid

    logical,parameter :: ifix_mass_col=.true.
!!    logical,parameter :: ifix_con_col=.false.
    logical,parameter :: ifix_con_col=.true.
!!!    logical,parameter :: ifix_vol_col=.false.
    logical,parameter :: ifix_vol_col=.true.

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! set up the maximum number of iteration
    NITER = 300
    NITER_TOTAL = 4
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Check the mass transfer
!!!    write(*,*) "bf repair debug 1"
    if(debug) then
      do n=1,ag%L
        rr_bef(n)=0.0
        cr_bef(n)=0.0
        cs_bef(n)=0.0
        uactr(n)=0.0_PS
        udepr(n)=0.0_PS
        garml(n)=0.0_PS
        lsmltl(n)=0.0_PS
        udepar(n)=0.0_PS
        uactar(n)=0.0_PS
        uimmar(n)=0.0_PS

        rs_bef(n)=0.0
        uacti(n)=0.0_PS
        udepi(n)=0.0_PS
        udhfi(n)=0.0_PS
        uimmai(n)=0.0_PS
        garm(n)=0.0_PS
        lsmlt(n)=0.0_PS
        udepai(n)=0.0_PS
      enddo

      do i=1,gr%N_BIN
        do n=1,gr%L
          rr_bef(n)=rr_bef(n)+gr%MS(i,n)%mass(rmt)-gr%MS(i,n)%mass(rmat)
          uactr(n)=uactr(n)+gr%MS(i,n)%dmassdt(rmt,3)-gr%MS(i,n)%dmassdt(rmat,3)
          udepr(n)=udepr(n)+gr%MS(i,n)%dmassdt(rmt,1)-gr%MS(i,n)%dmassdt(rmat,1)
          garml(n)=garml(n)+gr%MS(i,n)%dmassdt(rmt,4)-gr%MS(i,n)%dmassdt(rmat,4)
          lsmltl(n)=lsmltl(n)+gr%MS(i,n)%dmassdt(rmt,10)-gr%MS(i,n)%dmassdt(rmat,10)

          cr_bef(n)=cr_bef(n)+gr%MS(i,n)%con
          udepar(n)=udepar(n)+gr%MS(i,n)%dmassdt(rmat,1)
          uactar(n)=uactar(n)+gr%MS(i,n)%dmassdt(rmat,3)
          uimmar(n)=uimmar(n)+gr%MS(i,n)%dmassdt(rmat,11)
        enddo
      enddo

      do i=1,gs%N_BIN
        do n=1,gs%L
          rs_bef(n)=rs_bef(n)+gs%MS(i,n)%mass(imt)-gs%MS(i,n)%mass(imat)
          uacti(n)=uacti(n)+gs%MS(i,n)%dmassdt(imt,7)-gs%MS(i,n)%dmassdt(imat,7)
          udhfi(n)=udhfi(n)+gs%MS(i,n)%dmassdt(imt,3)-gs%MS(i,n)%dmassdt(imat,3)
          udepi(n)=udepi(n)+gs%MS(i,n)%dmassdt(imt,1)-gs%MS(i,n)%dmassdt(imat,1)
          garm(n)=garm(n)+gs%MS(i,n)%dmassdt(imt,4)-gs%MS(i,n)%dmassdt(imat,4)
          lsmlt(n)=lsmlt(n)+gs%MS(i,n)%dmassdt(imt,10)-gs%MS(i,n)%dmassdt(imat,10)
          udepai(n)=udepai(n)+gs%MS(i,n)%dmassdt(imat,1)
          cs_bef(n)=cs_bef(n)+gs%MS(i,n)%con
          uimmai(n)=uimmai(n)+gs%MS(i,n)%dmassdt(imat,11)
        enddo
      enddo

      do n=1,gs%L
        rr_bef(n)=rr_bef(n)/ag%TV(n)%den
        udepr(n)=udepr(n)*gr%dt
        uactr(n)=uactr(n)*gr%dt
        garml(n)=garml(n)*gr%dt
        lsmltl(n)=lsmltl(n)*gr%dt
        udepar(n)=udepar(n)*gr%dt
        uactar(n)=uactar(n)*gr%dt
        uimmar(n)=uimmar(n)*gr%dt

        rs_bef(n)=rs_bef(n)/ag%TV(n)%den
        udepi(n)=udepi(n)*gr%dt
        uacti(n)=uacti(n)*gr%dt
        uimmai(n)=uimmai(n)*gr%dt
        udhfi(n)=udhfi(n)*gr%dt
        garm(n)=garm(n)*gr%dt
        lsmlt(n)=lsmlt(n)*gr%dt
        udepai(n)=udepai(n)*gr%dt
      enddo

      do ica=1,ncat_a
        do n=1,ag%L
          udepa(ica,n)=0.0_PS
          uacta(ica,n)=0.0_PS
        enddo
      enddo
!      do n=1,LMAX*ncamx
!        udepa(n,1)=0.0_PS
!      enddo
      do i=1,ga(1)%N_BIN
        do ica=1,ncat_a
          do n=1,ag%L
            udepa(ica,n)=udepa(ica,n)+ga(ica)%MS(i,n)%dmassdt(amt,1)
            uacta(ica,n)=uacta(ica,n)+ga(ica)%MS(i,n)%dmassdt(amt,3)
          enddo
        enddo
      end do
!      do n=1,LMAX*ncamx
!        udepa(n,1)=udepa(n,1)*gr%dt
!      enddo
      do ica=1,ncat_a
        do n=1,ag%L
          udepa(ica,n)=udepa(ica,n)*gr%dt
          uacta(ica,n)=uacta(ica,n)*gr%dt
        enddo
      enddo

!      do in=1,gs%n_bin*gs%L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
      do n = 1, gs%L
      do i = 1, gs%N_BIN
!        ierror1(in)=0
         if((gs%MS(i,n)%dmassdt(imat,9)>1.0e+03.or.gs%MS(i,n)%dmassdt(imat,9)<-1.0e+03)&
        .or.(gs%MS(i,n)%dmassdt(imat,9)>0.0.and.gs%MS(i,n)%dmassdt(imat,9)<=0.0))then
            write(*,'("dmdt5 0",3I4,20ES15.6)') KD(n),ID(n),JD(n),gs%MS(i,n)%dmassdt(imat,9)
         endif
      enddo
      enddo

      do n=1,ag%L
!        ierror1(n)=0
         if( ( ag%TV(n)%s_v(1)>0.10.or.ag%TV(n)%s_v_n(1)>0.10 ) &
              .or. (uactar(n)>0.0.and.KD(n)>=90) ) then
            write(fid_alog,*) "repair called from:",trim(from),itcall
            write(fid_alog,*) "repair 0:k,i,j,udepr,uactr,udepi,uacti,udhfi,lsmlt,garm,sv1,sv2,dsv1,dsv2"
            write(fid_alog,*) &
                 KD(n),ID(n),JD(n),udepr(n),uactr(n),udepi(n),uacti(n),udhfi(n),lsmlt(n),garm(n),lsmltl(n),garml(n),&
                 ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
            write(fid_alog,*) "repair 02:k,i,j,udepa(cat),udepar,udepai"
            write(fid_alog,*) &
                 KD(n),ID(n),JD(n),udepa(1:ncat_a,n),udepar(n),udepai(n)
            write(fid_alog,*) "repair 03:k,i,j,cr,cs",cr_bef(n),cs_bef(n)
            write(fid_alog,*) "repair 03:k,i,j,dcondt1_s", &
                 gs%MS(1:gs%N_BIN,n)%dcondt(1)
            write(fid_alog,*) "repair 03:k,i,j,dcondt7_s", &
                 gs%MS(1:gs%N_BIN,n)%dcondt(7)
!            write(fid_alog,*) "repair 03:k,i,j,dcondt1", &
!                gr%MS(1:gr%N_BIN,n)%dcondt(1)
!            write(fid_alog,*) "repair 03:k,i,j,dcondt3", &
!                gr%MS(1:gr%N_BIN,n)%dcondt(3)
            write(fid_alog,*) "qv",KD(n),ID(n),JD(n),ag%TV(n)%rv,ag%TV(n)%t
            write(fid_alog,*) "repair 05:",uacta(1:ncat_a,n),uactar(n)
            write(fid_alog,*) "repair 06:",uimmar(n),uimmai(n)
         endif
!dbg        if(sum(gr%MS(:,n)%dmassdt(rmt,1))<0.0) then
!dbg          ierror1(n)=1
!dbg        endif
!dbg        if(sum(gs%MS(:,n)%dmassdt(imt,7))>0.0) then
!dbg          ierror1(n)=1
!dbg        endif
!dbg        if(sum(gr%MS(:,n)%dmassdt(rmt,3))>0.0) then
!dbg          ierror1(n)=1
!dbg        endif
      enddo

      do n=1,ag%L
!        ierror1(n)=0
         if(qtp(n)<rr_bef(n)+rs_bef(n)) then
            write(*,197) KD(n),ID(n),JD(n),rr_bef(n),rs_bef(n),rr_bef(n)+rs_bef(n),qtp(n)
197         format("something wrong bef: n,rr,rs,rr+rs,qtp at",3I5,4ES15.6)
         endif
      enddo

      do n=1,ag%L
!        ierror1(n)=0
         if(ga(2)%MS(1,n)%dcondt(1)>0.5_PS) then
            write(*,22) gs%IS(1,n)%sh_type,gs%IS(1,n)%habit,gs%MS(1,n)%dcondt(1),gs%MS(1,n)%dcondt(2)
21          format("ap2:n,dcondt2,sv1,sv2",I5,3ES15.6)
22          format("ap2:shtype,habit,dcondt_vp,col",2I5,2ES15.6)
         endif
      enddo

      do n=1,ag%L
        icond1(n)=0
        if(abs(ga(1)%MS(1,n)%dmassdt(amt,3))>1.0e-30_PS) then
          icond1(n)=1
        endif
      enddo
      do i=1,gr%n_bin
        do n=1,ag%L
          if(gr%MS(i,n)%dmassdt(rmt,3)>1.0e-30_PS.or.gr%MS(i,n)%dmassdt(rmat,3)>1.0e-30) then
            icond1(n)=icond1(n)+2
          endif
        enddo
      enddo
!      if(any(icond1(1:ag%L)>0)) then
!        do n=1,ag%L
!          if(icond1(n)>0) then
!            write(*,*) "repair 1, k,i,j:",KD(n),ID(n),JD(n)
!            write(*,*) "n,svw,svi",n,ag%TV(n)%s_v(1),ag%TV(n)%s_v(2)
!            write(*,*) ga(1)%MS(1,n)%dmassdt(amt,3),ga(1)%MS(1,n)%dmassdt(ams,3),ga(1)%MS(1,n)%dmassdt(ami,3)
!            write(*,*) ga(1)%MS(1,n)%mass
!            write(*,*) ga(1)%MS(1,n)%con
!          end if
!        enddo
!        stop
!      endif

   endif ! debug

!dbg    do in=1,gs%N_BIN*gs%L
!dbg      n=(in-1)/gs%N_BIN+1
!dbg      i=in-(n-1)*gs%N_BIN
!dbg      if(KD(n)==21.and.JD(n)==4) then
!dbg        write(*,*) "repair check0:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%dvoldt(iacr,1) &
!dbg                      ,gs%MS(i,n)%dvoldt(iccr,1)
!dbg
!dbg      endif
!dbg    enddo

!!!    write(*,*) "af repair debug 1"

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !
    ! Iteration process
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if(irep_vap) then

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! vapor adjustment
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! initialization

      do j=1,mxntend
!        do in=1,mxnbin*LMAX
!          n=(in-1)/mxnbin+1
!          i=in-(n-1)*mxnbin
         do n = 1, LMAX
         do i = 1, mxnbin
          acc_mod_r(i,j,n)=1.0_PS
          acc_mod_s(i,j,n)=1.0_PS
        enddo
        enddo
      enddo
      do ica=1,ncamx
        do j=1,mxntend
!          do in=1,mxnbina*LMAX
!            n=(in-1)/mxnbina+1
!            i=in-(n-1)*mxnbina
           do n = 1, LMAX
           do i = 1, mxnbina
            acc_mod_a(i,j,ica,n)=1.0_PS
           enddo
           enddo
        enddo
      enddo
      do n=1,LMAX
        mark_mass_r(n) = 0
        mark_mass_s(n) = 0
      enddo
!      do in=1,ncamx*LMAX
!        n=(in-1)/ncamx+1
!        ica=in-(n-1)*ncamx
      do n = 1, LMAX
      do ica = 1, ncamx
        mark_mass_a(ica,n)=0
      enddo
      enddo
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c          write(*,*) "front of mass_bug",n,n

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do itr1 = 1, NITER
        ! calculate the budget of each group and modify the mass tendency
        call cal_mass_budget_vapor( ag, gr, gs, ncat_a,ga, mes_rc, &
               acc_mod_r, acc_mod_s, acc_mod_a, &
               mark_mod, &
               mark_mass_r, mark_mass_s, mark_mass_a, flagp_r, flagp_s, flagp_a,&
               ID,JD,KD)

        if( any(mark_mod(1:ag%L)>0) )  then
!!c             if( i==1 ) then
!!c                write(*,*) " repair > vapor budget is fine at grid:",itr
!!c             else
!!c                write(*,*) " repair > required iteration for vapor budget is ", i
!!c             end if
        else
          exit
        end if
      end do
      if( any(mark_mod(1:ag%L) > 0) ) then
        do n=1,ag%L
          if(mark_mod(n)>0) then
             LOG_ERROR("repair",*) "Warning: The vapor mass tendency was not repaired successfully."
             LOG_ERROR_CONT(*) "         at",KD(n),ID(n),JD(n)
             LOG_ERROR_CONT(*) "sv1,rv,rv1",ag%TV(n)%s_v(1),ag%TV(n)%rv,ag%TV(n)%rv_sat(1)
             LOG_ERROR_CONT(*) "W,mark_mod,mark_acc_r,s,a",ag%TV(n)%W,mark_mod(n), &
               mark_mass_r(n),mark_mass_s(n),mark_mass_a(1:ncat_a,n)
             LOG_ERROR_CONT(*) "dmdt1",(gr%MS(i,n)%dmassdt(rmt,1),i=1,gr%N_BIN)
             LOG_ERROR_CONT(*) "dmdt3",(gr%MS(i,n)%dmassdt(rmt,3),i=1,gr%N_BIN)
             LOG_ERROR_CONT(*) "accmodr_1",(acc_mod_r(i,1,n),i=1,gr%N_BIN)
             LOG_ERROR_CONT(*) "accmodr_3",(acc_mod_r(i,3,n),i=1,gr%N_BIN)
             LOG_ERROR_CONT(*) "accmoda1",acc_mod_a(1,1,1,n),acc_mod_a(1,3,1,n)
             call PRC_abort
          endif
        enddo
      end if
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! modify other tendencies
      if( any(mark_mass_r(1:gr%L) /= 0 )) then
        call mod_other_tendency_vap(level, gr, acc_mod_r, 1)
!!c          write(*,*) " repair > tendencies of rain were repaired."
      end if
      if( any(mark_mass_s(1:gs%L) /= 0 )) then
        call mod_other_tendency_vap(level, gs, acc_mod_s, 1)
!!c          write(*,*) " repair > tendencies of solid hydro were repaired."
      end if
      if(level>=4) then
        do ica=1,ncat_a
          if( any(mark_mass_a(ica,1:ga(ica)%L) /= 0) ) then
            call mod_other_tendency_ap_vap( ga, acc_mod_a, ica,1)
!!c                write(*,*) " repair > tendencies of aerosols were repaired.",ica
          end if
        end do
      end if
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end if

!dbg    do in=1,gs%N_BIN*gs%L
!dbg      n=(in-1)/gs%N_BIN+1
!dbg      i=in-(n-1)*gs%N_BIN
!dbg      if(KD(n)==21.and.JD(n)==4) then
!dbg        write(*,*) "repair check1:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%dvoldt(iacr,1) &
!dbg                      ,gs%MS(i,n)%dvoldt(iccr,1)
!dbg
!dbg      endif
!dbg    enddo

!!!    write(*,*) "af repair vap"

    if(irep_col) then


      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Adjustment for collision (interacting) tendency
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do itr2 = 1, NITER_TOTAL

      if(ifix_mass_col) then
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! mass adjustment by collision tendency
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! initialization
        do j=1,mxntend
!          do in=1,mxnbin*LMAX
!            n=(in-1)/mxnbin+1
!            i=in-(n-1)*mxnbin
           do n = 1, LMAX
           do i = 1, mxnbin
            acc_mod_r(i,j,n)=1.0_PS
            acc_mod_s(i,j,n)=1.0_PS
          enddo
          enddo
        enddo
        do ica=1,ncamx
          do j=1,mxntend
!            do in=1,mxnbina*LMAX
!              n=(in-1)/mxnbina+1
!              i=in-(n-1)*mxnbina
             do n = 1, LMAX
             do i = 1, mxnbina
              acc_mod_a(i,j,ica,n)=1.0_PS
            enddo
            enddo
          enddo
        enddo
        do n=1,LMAX
          mark_mass_r(n) = 0
          mark_mass_s(n) = 0
        enddo
!        do in=1,ncamx*LMAX
!          n=(in-1)/ncamx+1
!          ica=in-(n-1)*ncamx
        do n = 1, LMAX
        do ica = 1, ncamx
          mark_mass_a(ica,n)=0
        enddo
        enddo
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! mass tendency
        do itr1 = 1, NITER
          ! calculate the budget of each group and modify the mass tendency
          call cal_mass_budget_col( ag, gr, gs, ncat_a,ga, &
                  acc_mod_r, acc_mod_s, acc_mod_a, &
                  mark_mod, &
                  mark_mass_r, mark_mass_s, mark_mass_a, flagp_r, flagp_s, flagp_a,&
                  iupdate_gr,iupdate_gs,iupdate_ga&
                  )
          if( any(mark_mod(1:ag%L) > 0 )) then

          else
!!c                if( i==1 ) then
!!c                   write(*,*) " repair > mass budget is fine."
!!c                else
!!c                   write(*,*) " repair > required iteration for mass budget is ", i
!!c                end if
            exit
          end if
        end do
!!!    write(*,*) "af repair col 1"
        if( any(mark_mod(1:ag%L)>0) )then
          do n=1,ag%L
            if(mark_mod(n)>0) then
               LOG_ERROR("repair",*) "Warning: The mass tendency was not repaired successfully."
               LOG_ERROR_CONT(*) "         at",KD(n),ID(n),JD(n)
               LOG_ERROR_CONT(*) "sv1,rv,rv1",ag%TV(n)%s_v(1),ag%TV(n)%rv,ag%TV(n)%rv_sat(1)
               LOG_ERROR_CONT(*) "W",ag%TV(n)%W
               LOG_ERROR_CONT(*) "dmdt1",(gr%MS(i,n)%dmassdt(rmt,1),i=1,gr%N_BIN)
               LOG_ERROR_CONT(*) "dmdt3",(gr%MS(i,n)%dmassdt(rmt,3),i=1,gr%N_BIN)
               LOG_ERROR_CONT(*) "accmodr_1",(acc_mod_r(i,1,n),i=1,gr%N_BIN)
               LOG_ERROR_CONT(*) "accmodr_3",(acc_mod_r(i,3,n),i=1,gr%N_BIN)
               LOG_ERROR_CONT(*) "accmoda1",acc_mod_a(1,1,1,n),acc_mod_a(1,3,1,n)
               call PRC_abort
            endif
          enddo
        end if
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! modify other tendencies
        if( any(mark_mass_r(1:gr%L) /= 0 )) then
          call mod_other_tendency(gr, acc_mod_r, 1, iupdate_gr)
!!c             write(*,*) " repair > tendencies of rain were repaired."
        end if
        if( any(mark_mass_s(1:gs%L) /= 0 )) then
          call mod_other_tendency(gs, acc_mod_s, 1, iupdate_gs)
!!c             write(*,*) " repair > tendencies of solid hydro were repaired."
        end if
        if(level>=4) then
          do ica=1,ncat_a
            if( any(mark_mass_a(ica,1:ga(ica)%L) /= 0) ) then
              call mod_other_tendency_ap( ga, acc_mod_a, ica,1, iupdate_ga)
!!c             write(*,*) " repair > tendencies of cloud droplets were repaired."
            end if
          end do
        end if
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      endif
!!!    write(*,*) "af repair col 2"

      if(ifix_con_col) then

        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! initialization
        do j=1,mxntend
!          do in=1,mxnbin*LMAX
!            n=(in-1)/mxnbin+1
!            i=in-(n-1)*mxnbin
           do n = 1, LMAX
           do i = 1, mxnbin
            acc_mod_r(i,j,n)=1.0_PS
            acc_mod_s(i,j,n)=1.0_PS
          enddo
          enddo
        enddo
        do ica=1,ncamx
          do j=1,mxntend
!            do in=1,mxnbina*LMAX
!              n=(in-1)/mxnbina+1
!              i=in-(n-1)*mxnbina
             do n = 1, LMAX
             do i = 1, mxnbina
              acc_mod_a(i,j,ica,n)=1.0_PS
            enddo
            enddo
          enddo
        enddo
        do n=1,LMAX
          mark_con_r(n) = 0
          mark_con_s(n) = 0
        enddo
!        do in=1,ncamx*LMAX
!          n=(in-1)/ncamx+1
!          ica=in-(n-1)*ncamx
        do n = 1, LMAX
        do ica = 1, ncamx
          mark_con_a(ica,n)=0
        enddo
        enddo
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! concentration tendency
        do itr1 = 1, NITER
          ! calculate the budget of each group and modify the mass tendency
          call cal_con_budget_col( gr, gs, ncat_a,ga, &
                  acc_mod_r, acc_mod_s,acc_mod_a, &
                  mark_mod, &
                  mark_con_r, mark_con_s, mark_con_a, flagp_r, flagp_s, flagp_a,&
                  iupdate_gr,iupdate_gs,iupdate_ga&
                  )
          if( any(mark_mod(1:ag%L) > 0 )) then

          else
!!c                if(i==1) then
!!c                   write(*,*) " repair > con budget is fine."
!!c                else
!!c                   write(*,*) " repair > required iteration for con budget is ", i
!!c                end if
            exit
          end if
        end do
!!!    write(*,*) "af repair col 3"
        if ( debug ) then
           if( any(mark_mod(1:ag%L)>0) )then
              do n=1,ag%L
                 if(mark_mod(n)>0) then
                    write(*,*) "Warning: The con tendency was not repaired successfully."
                    write(*,*) "         at",KD(n),ID(n),JD(n)
!!c             stop
                 end if
              enddo
           endif
        end if
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! modify other tendencies
        if( any(mark_con_r(1:gr%L) /= 0 )) then
          call mod_other_tendency(gr, acc_mod_r, 2, iupdate_gr)
!!c             write(*,*) " repair > tendencies of rain were repaired."
        end if
        if( any(mark_con_s(1:gs%L) /= 0 )) then
          call mod_other_tendency(gs, acc_mod_s, 2, iupdate_gs)
!!c             write(*,*) " repair > tendencies of solid hydro were repaired."
        end if
        if(level>=4) then
          do ica=1,ncat_a
            if( any(mark_con_a(ica,1:ga(ica)%L) /= 0) ) then
              call mod_other_tendency_ap( ga, acc_mod_a, ica,2, iupdate_ga)
!!c             write(*,*) " repair > tendencies of cloud droplets were repaired."
            end if
          end do
        end if
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!    write(*,*) "af repair col 4"
      endif

      if(ifix_vol_col) then
!org        do nvol=1,3
        do nvol=1,1
!test        do nvol=1,0
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! NOTE
          !   fixing PPVs are hard since we need to keep the relationship
          !   between them in a bin.
          !   So, currently only Vcs is considered.
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! initialization
          do j=1,mxntend
!            do in=1,mxnbin*LMAX
!              n=(in-1)/mxnbin+1
!              i=in-(n-1)*mxnbin
             do n = 1, LMAX
             do i = 1, mxnbin
              acc_mod_r(i,j,n)=1.0_PS
              acc_mod_s(i,j,n)=1.0_PS
            enddo
            enddo
          enddo
          do ica=1,ncamx
            do j=1,mxntend
!              do in=1,mxnbina*LMAX
!                n=(in-1)/mxnbina+1
!                i=in-(n-1)*mxnbina
               do n = 1, LMAX
               do i = 1, mxnbina
                acc_mod_a(i,j,ica,n)=1.0_PS
              enddo
              enddo
            enddo
          enddo
          do n=1,LMAX
            mark_v_r(n) = 0
            mark_v_s(n) = 0
          enddo
!          do in=1,ncamx*LMAX
!            n=(in-1)/ncamx+1
!            ica=in-(n-1)*ncamx
          do n = 1, LMAX
          do ica = 1, ncamx
            mark_v_a(ica,n)=0
          enddo
          enddo
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! circumscribing-volume tendency
          do itr1 = 1, NITER
            ! calculate the budget of each group and modify the Va tendency
            call cal_vol_budget_col( gr, gs, nvol, &
                     acc_mod_r, acc_mod_s, acc_mod_a, &
                     mark_mod, &
                     mark_v_r, mark_v_s, mark_v_a, flagp_s,&
                     iupdate_gs&
                     )
            if( any(mark_mod(1:ag%L) > 0 )) then

            else
!!c                   if(i==1)then
!!c                      write(*,*) " repair > volume budget (",nvol,") is fine."
!!c                   else
!!c                      write(*,*) " repair > required iteration for volume budget (",nvol,") is ", i
!!c                   end if
              exit
            end if
          end do
!!!    write(*,*) "af repair col 5"
          if ( debug ) then
             if( any(mark_mod(1:ag%L)>0) )then
                do n=1,ag%L
                   if(mark_mod(n)>0) then
                      write(*,*) "Warning: The volume tendency was not repaired successfully."
                      write(*,*) "         at",KD(n),ID(n),JD(n)
                   end if
                enddo
             endif
          end if
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! modify other tendencies
          if( any(mark_v_r(1:gr%L) /= 0 )) then
            call mod_other_tendency(gr, acc_mod_r, 2+nvol, iupdate_gr)
!!c                write(*,*) " repair > tendencies of rain were repaired."
          end if
          if( any(mark_v_s(1:gs%L) /= 0 )) then
            call mod_other_tendency(gs, acc_mod_s, 2+nvol, iupdate_gs)
!!c                write(*,*) " repair > tendencies of solid hydro were repaired."
          end if
!tmp             if(level>=4) then
!tmp                do ica=1,ncat_a
!tmp                   if( mark_v_a(ica) /= 0 ) then
!tmp                      call mod_other_tendency_ap( ga, n, acc_mod_a,ica,2+nvol, iupdate_ga)
!tmp!!c             write(*,*) " repair > tendencies of cloud droplets were repaired."
!tmp                   end if
!tmp                end do
!tmp             end if
             ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        end do

      endif

        mark_all=mark_mass_r+mark_mass_s+ &
                 mark_mass_a(1,:)+mark_mass_a(2,:) +&
                 mark_con_a(1,:)+mark_con_a(2,:)   +&
                 mark_con_r + mark_con_s +&
                 mark_v_a(1,:) + mark_v_a(2,:) + &
                 mark_v_r + mark_v_s
        if(sum(mark_all)==0) then
!tmp          write(*,*) " repair > required total iteration was ", itr2
          exit
        end if
      end do
    endif
!!!    write(*,*) "af repair col"

!dbg    do in=1,gs%N_BIN*gs%L
!dbg      n=(in-1)/gs%N_BIN+1
!dbg      i=in-(n-1)*gs%N_BIN
!dbg      if(KD(n)==12.and.JD(n)==77) then
!dbg        write(*,*) "repair check2:",i,KD(n),ID(n),JD(n),gs%MS(i,n)%dvoldt(iacr,4) &
!dbg                      ,gs%MS(i,n)%dvoldt(iccr,4)
!dbg
!dbg      endif
!dbg    enddo

    if(debug) then
      do n=1,ag%L
        rr(n)=0.0
        cr(n)=0.0
        uactr2(n)=0.0_PS
        udepr2(n)=0.0_PS
        garml2(n)=0.0_PS
        lsmltl2(n)=0.0_PS
        udepar2(n)=0.0_PS
        uactar2(n)=0.0_PS
        uimmar2(n)=0.0_PS

      enddo
      do i=1,gr%N_BIN
        do n=1,gr%L
          NDXDT(n) = 0.0_PS
          NDXDT2(n) = 0.0_PS
        enddo
        do j = 1, gr%N_tendpros
          do n=1,gr%L
!            ierror1(n)=0
            NDXDT(n) = NDXDT(n)+gr%MS(i,n)%dmassdt(rmt,j)
            NDXDT2(n) = NDXDT2(n)+gr%MS(i,n)%dmassdt(rmat,j)
            if(gr%MS(i,n)%dmassdt(rmt,j)>1.0e+03.or.&
               gr%MS(i,n)%dmassdt(rmt,j)<-1.0e+03)then
               LOG_ERROR("repair",201) n,i,j,gr%MS(i,n)%dmassdt(rmt,j)
201            format("mt_r>grid,bin,process,dmassdt from repair",3I5,ES15.6)
               LOG_ERROR_CONT(*) "after"

                do k=1,gs%N_BIN
                   LOG_ERROR_CONT(*) k,gs%MS(k,n)%dmassdt(imt,j)
                end do

                do k=1,gr%N_BIN
                   LOG_ERROR_CONT(*) k,gr%MS(k,n)%dmassdt(rmt,j)
                end do
                call PRC_abort
            endif
          enddo
        enddo
      enddo
      do i=1,gr%N_BIN
        do n=1,gr%L
          cr(n)=cr(n)+gr%MS(i,n)%con
          uactr2(n)=uactr2(n)+gr%MS(i,n)%dmassdt(rmt,3)-gr%MS(i,n)%dmassdt(rmat,3)
          udepr2(n)=udepr2(n)+gr%MS(i,n)%dmassdt(rmt,1)-gr%MS(i,n)%dmassdt(rmat,1)
          garml2(n)=garml2(n)+gr%MS(i,n)%dmassdt(rmt,4)-gr%MS(i,n)%dmassdt(rmat,4)
          lsmltl2(n)=lsmltl2(n)+gr%MS(i,n)%dmassdt(rmt,10)-gr%MS(i,n)%dmassdt(rmat,10)

          udepar2(n)=udepar2(n)+gr%MS(i,n)%dmassdt(rmat,1)
          uactar2(n)=uactar2(n)+gr%MS(i,n)%dmassdt(rmat,3)
          uimmar2(n)=uimmar2(n)+gr%MS(i,n)%dmassdt(rmat,11)

          nxr=max(0.0_PS,&
               (gr%MS(i,n)%mass(rmt)+NDXDT(n)*gr%dt-(gr%MS(i,n)%mass(rmat)+NDXDT2(n)*gr%dt))/ag%TV(n)%den)
!!c          if(gr%MS(i,n)%mass(rmt)+NDXDT*gr%dt<0.0_PS) then
!!c             write(*,*) "mass_r becomes negative",gr%MS(i,n)%mass(rmt)
!!c          end if

          rr(n)=rr(n)+nxr
        end do
      end do
      do n=1,gr%L
        uactr2(n)=uactr2(n)*gr%dt
        udepr2(n)=udepr2(n)*gr%dt
        garml2(n)=garml2(n)*gr%dt
        lsmltl2(n)=lsmltl2(n)*gr%dt
        udepar2(n)=udepar2(n)*gr%dt
        uactar2(n)=uactar2(n)*gr%dt
        uimmar2(n)=uimmar2(n)*gr%dt
      enddo

      if ( debug ) then
!      do in=1,gs%n_bin*gs%L
!        n=(in-1)/gs%N_BIN+1
!        i=in-(n-1)*gs%N_BIN
         do n = 1, gs%L
         do i = 1, gs%N_BIN
!        ierror1(in)=0
            if((gs%MS(i,n)%dmassdt(imat,9)>1.0e+03.or.gs%MS(i,n)%dmassdt(imat,9)<-1.0e+03)&
           .or.(gs%MS(i,n)%dmassdt(imat,9)>0.0.and.gs%MS(i,n)%dmassdt(imat,9)<=0.0))then
               write(*,'("dmdt5 1",3I4,20ES15.6)') KD(n),ID(n),JD(n),gs%MS(i,n)%dmassdt(imat,9)
            endif
         enddo
         enddo
      end if

      do n=1,ag%L
        rs(n)=0.0
        cs(n)=0.0
        uacti2(n)=0.0_PS
        udhfi2(n)=0.0_PS
        udepi2(n)=0.0_PS
        uimmai2(n)=0.0_PS
        garm2(n)=0.0_PS
        lsmlt2(n)=0.0_PS
        udepai2(n)=0.0_PS
      enddo
      do i=1,gs%N_BIN
        do n=1,gs%L
          NDXDT(n) = 0.0_PS
          NDXDT2(n) = 0.0_PS
        enddo
        do j = 1, gs%N_tendpros
          do n=1,gs%L
!            ierror1(n)=0
            NDXDT(n) = NDXDT(n)+gs%MS(i,n)%dmassdt(imt,j)
            NDXDT2(n) = NDXDT2(n)+gs%MS(i,n)%dmassdt(imat,j)
            if(gs%MS(i,n)%dmassdt(imt,j)>1.0e+03.or.&
               gs%MS(i,n)%dmassdt(imt,j)<-1.0e+03)then
               LOG_ERROR("repair",202) n,i,j,gs%MS(i,n)%dmassdt(imt,j)
202            format("mt_s>grid,bin,process,dmassdt from repair",3I5,ES15.6)
               LOG_ERROR_CONT(*) "after"
               LOG_ERROR_CONT(*) "ice"
               LOG_ERROR_CONT('(30ES15.6)') (gs%MS(k,n)%mass(1),k=1,gs%N_BIN)
                do k=1,gs%N_BIN
                   LOG_ERROR_CONT('(I5,15ES15.6)') k,(gs%MS(k,n)%dmassdt(imt,m),m=1,gs%N_tendpros)
                end do
!                LOG_ERROR_CONT(*) "bef tend of ice"
!                do k=1,gs%N_BIN
!                   LOG_ERROR_CONT('(I5,15ES15.6)') k,(temp_dmdt_s(k,n),n=1,gs%N_tendpros)
!                end do
                call PRC_abort
            endif
          enddo
        enddo
      enddo
      do i=1,gs%N_BIN
        do n=1,gs%L
          uacti2(n)=uacti2(n)+gs%MS(i,n)%dmassdt(imt,7)-gs%MS(i,n)%dmassdt(imat,7)
          udhfi2(n)=udhfi2(n)+gs%MS(i,n)%dmassdt(imt,3)-gs%MS(i,n)%dmassdt(imat,3)
          udepi2(n)=udepi2(n)+gs%MS(i,n)%dmassdt(imt,1)-gs%MS(i,n)%dmassdt(imat,1)
          garm2(n)=garm2(n)+gs%MS(i,n)%dmassdt(imt,4)-gs%MS(i,n)%dmassdt(imat,4)
          lsmlt2(n)=lsmlt2(n)+gs%MS(i,n)%dmassdt(imt,10)-gs%MS(i,n)%dmassdt(imat,10)
          udepai2(n)=udepai2(n)+gs%MS(i,n)%dmassdt(imat,1)
          uimmai2(n)=uimmai2(n)+gs%MS(i,n)%dmassdt(imat,11)

          nxs=max(0.0_PS,&
               (gs%MS(i,n)%mass(imt)+NDXDT(n)*gs%dt-(gs%MS(i,n)%mass(imat)+NDXDT2(n)*gs%dt))/ag%TV(n)%den)
          rs(n)=rs(n)+nxs

          cs(n)=cs(n)+gs%MS(i,n)%con
        end do
      end do
      do n=1,gs%L
        udepi2(n)=udepi2(n)*gr%dt
        uacti2(n)=uacti2(n)*gr%dt
        uimmai2(n)=uimmai2(n)*gr%dt
        udhfi2(n)=udhfi2(n)*gr%dt
        garm2(n)=garm2(n)*gr%dt
        lsmlt2(n)=lsmlt2(n)*gr%dt
        udepai2(n)=udepai2(n)*gr%dt
      enddo

!      do n=1,LMAX*ncamx
!        udepa2(n,1)=0.0_PS
!      enddo
      do ica=1,ncat_a
        do n=1,ag%L
          udepa2(ica,n)=0.0_PS
          uacta2(ica,n)=0.0_PS
        enddo
      enddo
      do i=1,ga(1)%N_BIN
        do ica=1,ncat_a
          do n=1,ag%L
            udepa2(ica,n)=udepa2(ica,n)+ga(ica)%MS(i,n)%dmassdt(amt,1)
            uacta2(ica,n)=uacta2(ica,n)+ga(ica)%MS(i,n)%dmassdt(amt,3)
          enddo
        enddo
      enddo
!      do n=1,LMAX*ncamx
!        udepa2(n,1)=udepa2(n,1)*gr%dt
!      enddo
      do ica=1,ncat_a
        do n=1,ag%L
          udepa2(ica,n)=udepa2(ica,n)*gr%dt
          uacta2(ica,n)=uacta2(ica,n)*gr%dt
        enddo
      enddo

!!c       write(*,*) "cr bf, cr",cr_bef,cr
!!c       write(*,'("bf dcdtr 1",25ES15.6)') (temp_dcdt_r(i,1),i=1,gr%N_BIN)
!!c       write(*,'("af dcdtr 1",25ES15.6)') (gr%MS(i,n)%dcondt(1),i=1,gr%N_BIN)
!!c       write(*,'("bf dcdtr 3",25ES15.6)') (temp_dcdt_r(i,3),i=1,gr%N_BIN)
!!c       write(*,'("af dcdtr 3",25ES15.6)') (gr%MS(i,n)%dcondt(3),i=1,gr%N_BIN)
!!c       write(*,'("bf dmdtr 1",25ES15.6)') (temp_dmdt_r(i,1),i=1,gr%N_BIN)
!!c       write(*,'("af dmdtr 1",25ES15.6)') (gr%MS(i,n)%dmassdt(rmt,1),i=1,gr%N_BIN)
!!c       write(*,'("bf dmdtr 3",25ES15.6)') (temp_dmdt_r(i,3),i=1,gr%N_BIN)
!!c       write(*,'("af dmdtr 3",25ES15.6)') (gr%MS(i,n)%dmassdt(rmt,3),i=1,gr%N_BIN)

!!c       write(*,198) n,rr,rs,rr+rs,qtp(n)
!!c198    format("aft: n,rr,rs,rr+rs,qtp at",I5,4ES15.6)
!!c       write(*,'(15ES15.6)') (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
!!c       write(*,'(15ES15.6)') (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)

      if ( debug ) then
         do n=1,ag%L
!        ierror1(n)=0
            if(qtp(n)<rr(n)+rs(n)) then
               write(*,198) n,rr(n),rs(n),rr(n)+rs(n),qtp(n)
198            format("something wrong aft: n,rr,rs,rr+rs,qtp at",I5,4ES15.6)
               write(*,*) "before rr,rs,rr+rs,qtp",rr_bef(n),rs_bef(n),rr_bef(n)+rs_bef(n),qtp(n)
               write(*,*) "density of air",ag%TV(n)%den
               write(*,*) "rain"
               write(*,'(15ES15.6)') (gr%MS(i,n)%mass(rmt),i=1,gr%N_BIN)
               do i=1,gr%N_BIN
                  write(*,'(I5,15ES15.6)') i,(gr%MS(i,n)%dmassdt(rmt,j),j=1,gr%N_tendpros)
               end do
!            write(*,*) "bef tend of rain"
!            do i=1,gr%N_BIN
!               write(*,'(I5,15ES15.6)') i,(temp_dmdt_r(i,j),j=1,gr%N_tendpros)
!            end do

               write(*,*) "ice"
               write(*,'(15ES15.6)') (gs%MS(i,n)%mass(imt),i=1,gs%N_BIN)
               do i=1,gs%N_BIN
                  write(*,'(I5,15ES15.6)') i,(gs%MS(i,n)%dmassdt(imt,j),j=1,gr%N_tendpros)
               end do
!            write(*,*) "bef tend of ice"
!            do i=1,gs%N_BIN
!               write(*,'(I5,15ES15.6)') i,(temp_dmdt_s(i,j),j=1,gs%N_tendpros)
!            end do
               write(*,'("sv1,sv2",5ES15.6)') ag%TV(n)%s_v(1),ag%TV(n)%s_v(2)
            endif
         enddo

         do n=1,ag%L
!        ierror1(n)=0
            if(ga(2)%MS(1,n)%dcondt(1)>0.5_PS) then
               write(*,11) n,ga(2)%MS(1,n)%dcondt(1),ag%TV(n)%s_v(1),ag%TV(n)%s_v(2)
               write(*,12) gs%IS(1,n)%sh_type,gs%IS(1,n)%habit,gs%MS(1,n)%dcondt(1),gs%MS(1,n)%dcondt(2)
11             format("ap3:n,dcondt2,sv1,sv2",I5,3ES15.6)
12             format("ap3:shtype,habit,dcondt_vp,col",2I5,2ES15.6)
            endif
         enddo
      end if ! debug

      do n=1,ag%L
        icond1(n)=0
        if(abs(ga(1)%MS(1,n)%dmassdt(amt,3))>1.0e-30_PS) then
          icond1(n)=1
        endif
      enddo
      do i=1,gr%n_bin
        do n=1,ag%L
          if(gr%MS(i,n)%dmassdt(rmt,3)>1.0e-30_PS.or.gr%MS(i,n)%dmassdt(rmat,3)>1.0e-30) then
            icond1(n)=icond1(n)+2
          end if
        end do
      end do
!      if(any(icond1(1:ag%L)==1)) then
      do n=1,ag%L
         if(icond1(n)==1) then

            LOG_ERROR("repair",*) "repair 2"
            LOG_ERROR_CONT(*) "n,svw,svi",n,ag%TV(n)%s_v(1),ag%TV(n)%s_v(2)
            LOG_ERROR_CONT(*) ga(1)%MS(1,n)%dmassdt(amt,3),ga(1)%MS(1,n)%dmassdt(ams,3),ga(1)%MS(1,n)%dmassdt(ami,3)
            LOG_ERROR_CONT(*) ga(1)%MS(1,n)%mass(1:3)
            LOG_ERROR_CONT(*) ga(1)%MS(1,n)%con
            call PRC_abort
         end if
      enddo
!      endif

      if ( debug ) then
         do n=1,ag%L
!        ierror1(n)=0
            if( (ag%TV(n)%s_v(1)>0.10.or.ag%TV(n)%s_v_n(1)>0.10) &
                 .or. (uactar(n)>0.0.and.KD(n)>=90) ) then
               write(fid_alog,*) "repair called from:",trim(from),itcall
               write(fid_alog,*) "repair 1:k,i,j,udepr,uactr,udepi,uacti,udhfi,lsmlt,garm,lsmltl,garml,sv1,sv2,dsv1,dsv2"
               write(fid_alog,*) &
                    KD(n),ID(n),JD(n),udepr2(n),uactr2(n),udepi2(n),uacti2(n),udhfi2(n),&
                    lsmlt2(n),garm2(n),lsmltl2(n),garml2(n),&
                    ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
               write(fid_alog,*) "repair 12:k,i,j,udepa(cat),udepar,udepai"
               write(fid_alog,*) &
                    KD(n),ID(n),JD(n),udepa2(1:ncat_a,n),udepar2(n),udepai2(n)
               write(fid_alog,*) "repair 13:k,i,j,cr,cs",cr(n),cs(n)
!            write(fid_alog),*) "repair 13:k,i,j,dcondt1", &
!                gr%MS(1:gr%N_BIN,n)%dcondt(1)
!            write(fid_alog),*) "repair 13:k,i,j,dcondt3", &
!                gr%MS(1:gr%N_BIN,n)%dcondt(3)
               write(fid_alog,*) "repair 13:k,i,j,dcondt1_s", &
                    gs%MS(1:gs%N_BIN,n)%dcondt(1)
               write(fid_alog,*) "repair 13:k,i,j,dcondt7_s", &
                    gs%MS(1:gs%N_BIN,n)%dcondt(7)
               write(fid_alog,*) "qv",KD(n),ID(n),JD(n),ag%TV(n)%rv,ag%TV(n)%t
               write(fid_alog,*) "repair 15:",uacta2(1:ncat_a,n),uactar2(n)
               write(fid_alog,*) "repair 16:",uimmar2(n),uimmai2(n)
            endif
!dbg        if(uacti(n)>0.0) then
!dbg          ierror1(n)=1
!dbg        endif
!dbg        if(udhfi(n)>0.0) then
!dbg          ierror1(n)=1
!dbg        endif
!dbg        if(uactr(n)>0.0) then
!dbg          ierror1(n)=1
!dbg        endif
!dbg        if(udepar(n)<0.0) then
!dbg          ierror1(n)=1
!dbg        endif
         enddo
      endif

   end if

!!!    write(*,*) "af repair debug 2"

!!c       write(*,'("repair:k,i,j,udepi,uacti,udepi2,uacti2,sv1,sv2,dsv1,dsv2",3I5,10ES15.6)') KD(n),ID(n),JD(n),udepi,uacti,udepi2,uacti2,ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)

!!c       if(abs(lsmlt2)>0.0.and.(abs(lsmltl2)-abs(lsmlt2))/abs(lsmlt2)>1.0e-4) then
!!c          write(*,*) "difference in water transfer between liq and ice"
!!c          write(*,'("repair 0:k,i,j,udepr,uactr,udepi,uacti,lsmlt,garm,sv1,sv2,dsv1,dsv2",3I5,30ES15.6)') &
!!c            KD(n),ID(n),JD(n),udepr,uactr,udepi,uacti,lsmlt,garm,lsmltl,garml,&
!!c            ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
!!c
!!c          write(*,'("repair:k,i,j,udepr,uactr,udepi,uacti,lsmlt,garm,lsmltl,garml,sv1,sv2,dsv1,dsv2",3I5,30ES15.6)') &
!!c          KD(n),ID(n),JD(n),udepr2,uactr2,udepi2,uacti2,lsmlt2,garm2,lsmltl2,garml2,&
!!c          ag%TV(n)%s_v(1),ag%TV(n)%s_v(2),ag%TV(n)%s_v_n(1),ag%TV(n)%s_v_n(2)
!!c
!!c                write(*,*) "af tend of ice"
!!c                do k=1,gs%N_BIN
!!c                   write(*,'(I5,15ES15.6)') k,(gs%MS(k,n)%dmassdt(imt,n),n=1,gs%N_tendpros)
!!c                end do
!!c                write(*,*) "bf tend of ice"
!!c                do k=1,gs%N_BIN
!!c                   write(*,'(I5,15ES15.6)') k,(temp_dmdt_s(k,n),n=1,gs%N_tendpros)
!!c                end do
!!c       endif

!!c       end if
!!c       write(*,22) n,rr,rs,rr+rs,qtp(n)
!!c       temp_dmdtap_s1=0.0_PS
!!c       do i=1,gs%N_BIN
!!c          temp_dmdtap_s1=temp_dmdtap_s1+gs%MS(i,n)%dmassdt(imat,1)
!!c       end do
!!c       if(temp_dmdtap_s1/=temp_dmdtap_s0) then
!!c          write(*,'("n,dmdtap0,dmdtap1",I5,2ES15.6)') n,temp_dmdtap_s0,temp_dmdtap_s1
!!c          write(*,*) "org tend"
!!c          do i=1,gs%N_BIN
!!c             write(*,'(20ES15.6)') (temp_dmdtap_s(i,j),j=1,gs%N_tendpros)
!!c          end do
!!c          write(*,*) "mod tend"
!!c          do i=1,gs%N_BIN
!!c             write(*,'(20ES15.6)') (gs%MS(i,n)%dmassdt(imat,j),j=1,gs%N_tendpros)
!!c          end do
!!c       end if
!!c       if(n==810.or.n==885.or.n==886) then
!!c          do i=1,gs%N_BIN
!!c          write(*,*) "rep 2:i,n",i,n
!!c          write(*,'("dmdt1",20ES15.6)') (gs%MS(i,n)%dmassdt(imt,j),j=1,gs%n_tendpros)
!!c          write(*,'("dcdt1",20ES15.6)') (gs%MS(i,n)%dcondt(j),j=1,gs%n_tendpros)
!!c          write(*,'("dvoldt1",20ES15.6)') (gs%MS(i,n)%dvoldt(ivcs,j),j=1,gs%n_tendpros)
!!c          write(*,'("dvoldt2",20ES15.6)') (gs%MS(i,n)%dvoldt(iacr,j),j=1,gs%n_tendpros)
!!c          write(*,'("dvoldt3",20ES15.6)') (gs%MS(i,n)%dvoldt(iccr,j),j=1,gs%n_tendpros)
!!c          end do
!!c       end if
!tmp    call deallocate_repair
!tmp  contains
!tmp    subroutine allocate_repair
!tmp      integer :: item(6)
!tmp      item=0
!tmp      allocate( acc_mod_r(gr%N_BIN, gr%N_tendpros), stat = item(1))
!tmp      allocate( acc_mod_s(gs%N_BIN, gs%N_tendpros), stat = item(2))
!tmp      allocate( acc_mod_a(max(ga(1)%N_BIN,ga(2)%N_BIN),ga(1)%N_tendpros,2), stat = item(3))
!tmp      allocate( mark_mass_a(2), stat = item(4))
!tmp      allocate( mark_con_a(2), stat = item(5))
!tmp      allocate( mark_v_a(2), stat = item(6))
!tmp      if(any(item /= 0 ) ) then
!tmp         write(*,*) item
!tmp         stop "Memory not available in repair"
!tmp      end if
!tmp    end subroutine allocate_repair
!tmp    subroutine deallocate_repair
!tmp      integer :: item(6)
!tmp      item=0
!tmp      deallocate( mark_v_a, stat = item(6))
!tmp      nullify(mark_v_a)
!tmp      deallocate( mark_con_a, stat = item(5))
!tmp      nullify(mark_con_a)
!tmp      deallocate( mark_mass_a, stat = item(4))
!tmp      nullify(mark_mass_a)
!tmp      deallocate( acc_mod_a, stat = item(3))
!tmp      nullify(acc_mod_a)
!tmp      deallocate( acc_mod_s, stat = item(2))
!tmp      nullify(acc_mod_s)
!tmp      deallocate( acc_mod_r, stat = item(1))
!tmp      nullify(acc_mod_r)
!tmp      if(any(item /= 0 )) then
!tmp         write(*,*) item
!tmp         stop "Memory cannot be deallocated in repair"
!tmp      end if
!tmp    end subroutine deallocate_repair
  end subroutine repair
  

  ! private procedures

  subroutine cal_mass_budget_col( ag, gr, gs, ncat_a,ga, &
       acc_mod_r, acc_mod_s, acc_mod_a,&
       mark_mod, mark_acc_r , mark_acc_s, mark_acc_a, flagp_r, flagp_s, flagp_a,&
       iupdate_gr,iupdate_gs,iupdate_ga&
       )
    use scale_prc, only: &
       PRC_abort
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate total gain and loss for each group, and modify
    ! only mass tendency
    !
    ! Collision only
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (AirGroup), intent(inout)    :: ag
    type (Group), intent(inout)         :: gr, gs
    integer,intent(in) :: ncat_a
    type (Group), dimension(*)          :: ga

    ! message from reality_check subroutine
    ! 0. vapor and hydrometeors do not exist.
    ! 1. vapor only. no hydrometeors.
    ! 2. vapor and liquid hydrometeors exist.
    ! 3. vapor and solid hydrometeoors exist.
    ! 4. vapor, liquid and solid hydrometeors exist.
    !
    !integer,dimension(*)   :: mes_rc

    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_r
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_s
    real (PS), dimension(mxnbina,mxntend,ncamx,*)    :: acc_mod_a

    integer, dimension(*), intent(inout)     :: mark_mod
    integer, dimension(ncamx,*)    :: mark_acc_a ! mark_mass_a
    integer, dimension(*),intent(inout)     :: mark_acc_r, mark_acc_s ! mark_mass_r, mark_mass_s
    integer, intent(in)                   :: flagp_r, flagp_s, flagp_a
    !integer :: ID(*),JD(*),KD(*)
    integer,dimension(*),intent(in) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,*),intent(in) :: iupdate_ga
    !
    ! local space
    !
    ! mass of vapor (g/cm^3), saturation mass of vapor
    !real(PS),dimension(LMAX)   :: mass_v, mass_vs, mass_vc
    real(PS),dimension(LMAX)   :: total_rimr_tend,total_rims_tend,total_mlts_tend,total_mltr_tend
    real(PS)    :: total_rloss, total_rgain, total_sloss, total_sgain,total_aloss, total_again

    ! modification constant for vapor, cloud_drop, rain, and solid_hydro
    real(PS),dimension(LMAX)    :: modc_r, modc_s, modc_a!, modc_v
    real(PS),dimension(LMAX)    :: mod_dum,dum1,dum2


    ! temperature of freezing
    real(PS),parameter :: total_limit=1.0e-30_PS


    integer,dimension(LMAX) :: icond1,icond2,icond3 !,ierror1
    real(PS) :: exist_mass

    integer      :: i,j, k, ica,n

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize
    do n=1,ag%L
      mark_mod(n) = 0
    enddo
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2. budget of collision processes and others
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2.1. bugdet of each rain bin
    ! calculate total tendency of riming of rain drops by ices
    if( flagp_r > 0 ) then

      liqbin_loop1: do i=1,gr%n_bin
        do n=1,gr%L

          total_rloss = &
               ! - collision advection
               max(-gr%MS(i,n)%dmassdt(rmt,2),0.0_ds)*real(iupdate_gr(2),PS_KIND) + &
               ! - riming process by a solid hyrometeor
               max(-gr%MS(i,n)%dmassdt(rmt,4), 0.0_ds)*real(iupdate_gr(4),PS_KIND) + &
               ! - collision breakup
               max(-gr%MS(i,n)%dmassdt(rmt,5),0.0_ds)*real(iupdate_gr(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(-gr%MS(i,n)%dmassdt(rmt,6),0.0_ds)*real(iupdate_gr(6),PS_KIND) + &
!!c               ! - evaporation - moist adiabatic ascent - for cloud droplet bin
!!c               max(-gr%MS(i,n)%dmassdt(rmt,7),0.0_ds) + &
               ! - loss of cloud droplets by contact nucleation
               max(-gr%MS(i,n)%dmassdt(rmt,8),0.0_ds)*real(iupdate_gr(8),PS_KIND) + &
               ! - auto conversion
               max(-gr%MS(i,n)%dmassdt(rmt,9), 0.0_ds)*real(iupdate_gr(9),PS_KIND) + &
               ! - ice nucleation (immersion)
               max(-gr%MS(i,n)%dmassdt(rmt,11), 0.0_ds)*real(iupdate_gr(11),PS_KIND) + &
               ! - ice nucleation (homogeneous freezing)
               max(-gr%MS(i,n)%dmassdt(rmt,12), 0.0_ds)*real(iupdate_gr(12),PS_KIND)

          total_rgain = &
               ! - existing mass
               gr%MS(i,n)%mass(rmt)/gr%dt + &
!!c               0.9999*gr%MS(i,n)%mass(rmt)/gr%dt + &
               ! - collision advection
!!c               min((rr-gr%MS(i,n)%mass(rmt))/gr%dt,max(gr%MS(i,n)%dmassdt(rmt,2),0.0_ps)) + &
               max(gr%MS(i,n)%dmassdt(rmt,2),0.0_ds)*real(iupdate_gr(2),PS_KIND) + &
               ! - collision breakup
               max(gr%MS(i,n)%dmassdt(rmt,5),0.0_ds)*real(iupdate_gr(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(gr%MS(i,n)%dmassdt(rmt,6),0.0_ds)*real(iupdate_gr(6),PS_KIND) + &
!!c               ! - vapor deposition - moist adiabatic ascent - for cloud droplet bin
!!c               max(gr%MS(i,n)%dmassdt(rmt,7),0.0_ds) + &
               ! - auto conversion
               max(gr%MS(i,n)%dmassdt(rmt,9), 0.0_ds)*real(iupdate_gr(9),PS_KIND) + &
               ! - melting-shedding
               max(gr%MS(i,n)%dmassdt(rmt,10),0.0_ds)*real(iupdate_gr(10),PS_KIND)

!          ierror1(n)=0
          icond1(n)=0
          modc_r(n)=1.0_PS

          if( total_rgain < 0.99999_ps*total_rloss.and.total_rloss>total_limit) then

            icond1(n)=1

            modc_r(n) = min(1.0_ps,&
                      (total_rgain)/max(total_rloss,total_limit))

            modc_r(n)=max(0.0_PS,modc_r(n))

            mark_mod(n) = mark_mod(n) + 1
            mark_acc_r(n) = mark_acc_r(n) + 1

          else if( total_rgain == 0.0_ps .and. total_rloss == 0.0_ps ) then
            ! this mark should help having bogus tendency in
            ! calculating new tendency (cal_model_tendency)
            gr%MS(i,n)%mark = 4
          end if
          if(modc_r(n)>1.0e+4) then
             LOG_ERROR("cal_mass_budget_col",*) "something wrong",n,gr%MS(i,n)%dmassdt(rmt,1:gr%n_tendpros)
             call prc_abort
          endif
        enddo

        if(all(icond1(1:gr%L)==0)) cycle

        do j=1,gr%n_tendpros
          if(iupdate_gr(j)==0) cycle

          do n=1,gr%L
            icond2(n)=0
            if( gr%MS(i,n)%dmassdt(rmt,j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo

          if(j == 2 .or. j == 5 .or. j == 6.or.j==9) then
            ! collision-coalescence
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if(icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
          elseif(j==4.and.flagp_s>0) then
            ! riming process
            ! this decreases total transfer
            do n=1,gr%L
              if( icond2(n)==1) then
                gr%MS(i,n)%dmassdt(rmt,j)=gr%MS(i,n)%dmassdt(rmt,j)*modc_r(n)
                acc_mod_r(i,j,n)=acc_mod_r(i,j,n)*modc_r(n)
              endif
              total_rimr_tend(n)=0.0_ps
              total_rims_tend(n)=0.0_ps
            enddo
            do k=1,gr%N_bin
              do n=1,gr%L
                total_rimr_tend(n)=total_rimr_tend(n)-gr%MS(k,n)%dmassdt(rmt,j)
              enddo
            enddo
            do k=1,gs%N_bin
              do n=1,gs%L
                total_rims_tend(n)=total_rims_tend(n)+gs%MS(k,n)%dmassdt(imt,j)
              enddo
            enddo
            do n=1,gr%L
              mod_dum(n)=1.0_PS
              if( icond2(n)==1) then
                mod_dum(n)=max(0.0_PS,min(1.0_PS,total_rimr_tend(n)/max(total_limit,total_rims_tend(n))))
!tmp              ag%tv(n)%dmassdt_v(3)=ag%tv(n)%dmassdt_v(3)*modc_r(n)
              endif
              mark_acc_s(n)=mark_acc_s(n)+icond1(n)*icond2(n)
            enddo

!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
            do n = 1, gs%L
            do k = 1, gs%n_bin
              gs%MS(k,n)%dmassdt(imt,j)=gs%MS(k,n)%dmassdt(imt,j)*mod_dum(n)
              acc_mod_s(k,j,n)=acc_mod_s(k,j,n)*mod_dum(n)
            end do
            end do

          else if((j == 11 .or. j == 12 ).and. flagp_s > 0 ) then
            ! - immersion freezing or homogeneous freezing
            !   all the tendencies in those bins are reduced.
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
            do n = 1, gs%L
            do k = 1, gs%n_bin
              if( icond2(n)==1 ) then
                gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_r(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n=1,gs%L
              mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
            enddo
          else if( j == 8 ) then
            ! contact ice nucleation
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
            do n = 1, gs%L
            do k = 1, gs%n_bin
              if( icond2(n)==1) then
                gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_r(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n=1,gs%L
              mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
            enddo

            do n=1,ag%L
              if( icond2(n)==1) then
                ga(2)%ms(1,n)%dmassdt(amt,j) = ga(2)%ms(1,n)%dmassdt(amt,j) * modc_r(n)
                acc_mod_a(1,j,2,n) = acc_mod_a(1,j,2,n) * modc_r(n)
              endif
              mark_acc_a(2,n)=mark_acc_a(2,n) + icond1(n)*icond2(n)
            enddo
          else
            do n=1,ag%L
              if( icond2(n)==1) then
                gr%MS(i,n)%dmassdt(rmt,j) = gr%MS(i,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(i,j,n) = acc_mod_r(i,j,n) * modc_r(n)
              end if
            enddo
          endif
        end do
      enddo liqbin_loop1
    end if
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2.2. bugdet of each solid hydrometeor bin
    if( flagp_s > 0 ) then
      rimmlt_if: if(iupdate_gs(4)==1.and.iupdate_gs(10)==1) then
        ! first adjust melting mass
        do n=1,gr%L
          total_mltr_tend(n)=0.0_PS
          total_rimr_tend(n)=0.0_PS
          total_mlts_tend(n)=0.0_PS
          total_rims_tend(n)=0.0_PS
        enddo
        do k=1,gr%n_bin
          do n=1,gr%L
            total_mltr_tend(n)=total_mltr_tend(n)+gr%MS(k,n)%dmassdt(rmt,10)
            total_rimr_tend(n)=total_rimr_tend(n)-gr%MS(k,n)%dmassdt(rmt,4)
          enddo
        enddo
        do k=1,gs%n_bin
          do n=1,gs%L
            total_mlts_tend(n)=total_mlts_tend(n)-gs%MS(k,n)%dmassdt(imt,10)
            total_rims_tend(n)=total_rims_tend(n)+gs%MS(k,n)%dmassdt(imt,4)
          enddo
        enddo
!tmp        do n=1,gs%L
!tmp          modc_s(n)=total_mltr_tend(n)/max(1.0e-30_PS,total_mlts_tend(n))
!tmp        enddo

!!c       if(total_mlts_tend>0.0_PS.or.total_rims_tend>0.0_PS) then
!!c          write(*,'("ice>mltck1",3I5,20ES15.6)') KD(n),ID(n),JD(n),&
!!c                      modc_s,total_mltr_tend*gr%dt,total_mlts_tend*gr%dt,&
!!c                      total_rimr_tend*gr%dt,total_rims_tend*gr%dt
!!c       endif

        do n=1,gs%L
          modc_s(n)=1.0_PS
!tmp          if(total_mltr_tend(n)>total_limit.and.abs(modc_s(n)-1.0)>1.0e-5) then
          if(total_mltr_tend(n)<0.99999_PS*total_mlts_tend(n).and. &
             total_mlts_tend(n)>total_limit) then

            modc_s(n)=min(1.0_PS,&
                    total_mltr_tend(n)/max(total_mlts_tend(n),total_limit))

            mark_mod(n) = mark_mod(n) + 1
            mark_acc_s(n) = mark_acc_s(n) + 1
          endif
        enddo
!        do kn=1,gs%n_bin*gs%L
!          n=(kn-1)/gs%n_bin+1
!          k=kn-(n-1)*gs%n_bin
        do n = 1, gs%L
        do k = 1, gs%n_bin
          gs%MS(k,n)%dmassdt(imt,10)=gs%MS(k,n)%dmassdt(imt,10)*modc_s(n)
          acc_mod_s(k,10,n)=acc_mod_s(k,10,n)*modc_s(n)
        end do
        end do

        icebin_loop1: do i=1,gs%n_bin
          do n=1,gs%L
            total_sloss = &
               ! - riming process with rain drops
               max(-gs%MS(i,n)%dmassdt(imt,4),0.0_ds) + &
               ! - melting-shedding
               max(-gs%MS(i,n)%dmassdt(imt,10),0.0_ds)


            total_sgain = &
               ! - existing mass
               gs%MS(i,n)%mass(imt)/gs%dt + &
!!c               0.9999*gs%MS(i,n)%mass(imt)/gs%dt+&
               ! - riming process with rain drops
               max(gs%MS(i,n)%dmassdt(imt,4),0.0_ds)

            icond1(n)=0
!            ierror1(n)=0
            modc_s(n)=1.0_PS

            if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then

              modc_s(n) = min(1.0_ps,&
                total_sgain/max(total_sloss,total_limit))

!!!             if(total_sloss<total_limit) cycle

              modc_s(n)=max(0.0_PS,modc_s(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_s(n) = mark_acc_s(n) + 1

              icond1(n)=1

!!c             write(*,23) i,KD(n),ID(n),JD(n),modc_s,total_sgain,total_sloss
!!c23           format("mass:modc_s,gain,loss",4i5,3es15.6)
            endif
            if(modc_s(n)>1.0e+4) then
               LOG_ERROR("cal_mass_budget_col",*) "something wrong",gs%MS(i,n)%dmassdt(imt,1:gs%n_tendpros)
               call PRC_abort
            endif
          enddo

          if(all(icond1(1:gs%L)==0)) cycle

          do j=4,10,6
            do n=1,gs%L
              icond2(n)=0
              if(gs%MS(i,n)%dmassdt(imt,j) < 0.0_ps) then
                icond2(n)=1
              endif
            enddo

            if(j==4) then
              ! riming process with rain drops
              ! this tries to shrink positive tendency, so that total transfer
              ! from liq to ice does not change.
              !
              ! calculate total transfer from liq to ice
              do n=1,gs%L
                dum1(n)=0.0_ps
                dum2(n)=0.0_ps
              enddo

              do k=1,gs%n_bin
                do n=1,gs%L
                  if(gs%MS(k,n)%dmassdt(imt,j)>0.0_PS) then
                    dum1(n)=dum1(n)+gs%MS(k,n)%dmassdt(imt,j)*real(icond2(n),PS_KIND)
                  else
                    dum2(n)=dum2(n)-gs%MS(k,n)%dmassdt(imt,j)*real(icond2(n),PS_KIND)
                  endif
                enddo
              enddo
              do n=1,gs%L
                total_rims_tend(n)=dum1(n)-dum2(n)

                dum2(n)=dum2(n)+(modc_s(n)-1.0_PS)*(-gs%MS(i,n)%dmassdt(imt,j))*real(icond2(n),PS_KIND)

                mod_dum(n)=1.0_PS
                icond3(n)=0
                if(real(icond2(n),PS_KIND)*dum1(n)>1.0e-30_PS) then
                  mod_dum(n)=max(0.0_DS,min(1.0_DS,&
                       (dum1(n)+(modc_s(n)-1.0_PS)*(-gs%MS(i,n)%dmassdt(imt,j)))/dum1(n)))

                  gs%MS(i,n)%dmassdt(imt,j)=gs%MS(i,n)%dmassdt(imt,j)*modc_s(n)
                  acc_mod_s(i,j,n)=acc_mod_s(i,j,n)*modc_s(n)
                  icond3(n)=1
                endif
              enddo

!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond3(n)==1.and.gs%MS(k,n)%dmassdt(imt,j)>0.0_PS) then
                  gs%MS(k,n)%dmassdt(imt,j)=gs%MS(k,n)%dmassdt(imt,j)*mod_dum(n)
                  acc_mod_s(k,j,n)=acc_mod_s(k,j,n)*mod_dum(n)
                endif
                if(icond3(n)==0.and.icond2(n)==1) then
                  gs%MS(k,n)%dmassdt(imt,j)=gs%MS(k,n)%dmassdt(imt,j)*modc_s(n)
                  acc_mod_s(k,j,n)=acc_mod_s(k,j,n)*modc_s(n)
                endif
              enddo
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond3(n)==0.and.icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j)=gr%MS(k,n)%dmassdt(rmt,j)*modc_s(n)
                  acc_mod_r(k,j,n)=acc_mod_r(k,j,n)*modc_s(n)
                endif
              enddo
              enddo
              do n=1,gr%L
                mark_acc_r(n)=mark_acc_r(n)+icond1(n)*(1-icond3(n))*icond2(n)
              enddo

!!c                      total_rims_tend=0.0_ps
!!c                      do k=1,gs%n_bin
!!c                         total_rims_tend=total_rims_tend+gs%MS(k,n)%dmassdt(imt,4)
!!c                      end do
!!c                      total_rimr_tend=0.0_ps
!!c                      do k=1,gr%n_bin
!!c                         total_rimr_tend=total_rimr_tend-gr%MS(k,n)%dmassdt(rmt,4)
!!c                      end do
!!c                      write(*,'("ice>mod_dum,rim",4I5,10ES15.6)') KD(n),ID(n),JD(n),i,mod_dum,modc_s,&
!!c                                   total_rims_tend*gs%dt,total_rimr_tend*gs%dt

            elseif(j==10) then
              ! melting-shedding process
              ! total transfer from ice to liq does decrease.
              !
              ! calculate total transfer from ice to liq or liq to ice
              ! all tendencies for melting process is positive for liq and neg for ice.
              do n=1,gs%L
                if(icond2(n)==1) then
                  gs%MS(i,n)%dmassdt(imt,j)=gs%MS(i,n)%dmassdt(imt,j)*modc_s(n)
                  acc_mod_s(i,j,n)=acc_mod_s(i,j,n)*modc_s(n)
                endif
              enddo
              do n=1,gs%L
                total_mlts_tend(n)=0.0_PS
                total_mltr_tend(n)=0.0_PS
              enddo

              do k=1,gs%n_bin
                do n=1,gs%L
                  total_mlts_tend(n)=total_mlts_tend(n)-gs%MS(k,n)%dmassdt(imt,j)
                enddo
              enddo
              do k=1,gr%n_bin
                do n=1,gs%L
                  total_mltr_tend(n)=total_mltr_tend(n)+gr%MS(k,n)%dmassdt(rmt,j)
                enddo
              enddo
              do n=1,gs%L
                mod_dum(n)=max(0.0_PS,min(1.0_PS,total_mlts_tend(n)/ &
                           max(total_limit,total_mltr_tend(n))))
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j)=gr%MS(k,n)%dmassdt(rmt,j)*mod_dum(n)
                  acc_mod_r(k,j,n)=acc_mod_r(k,j,n)*mod_dum(n)
                endif
              enddo
              enddo
              do n=1,gs%L
                mark_acc_r(n) = mark_acc_r(n) + icond1(n)*icond2(n)
              enddo

            endif
          enddo

        enddo icebin_loop1
      endif rimmlt_if
!!c       total_mlts_tend=0.0_ps
!!c       total_rims_tend=0.0_ps
!!c       do k=1,gs%n_bin
!!c          total_mlts_tend=total_mlts_tend-gs%MS(k,n)%dmassdt(imt,10)
!!c          total_rims_tend=total_rims_tend+gs%MS(k,n)%dmassdt(imt,4)
!!c       end do
!!c       total_mltr_tend=0.0
!!c       total_rimr_tend=0.0_ps
!!c       do k=1,gr%n_bin
!!c          total_mltr_tend=total_mltr_tend+gr%MS(k,n)%dmassdt(rmt,10)
!!c          total_rimr_tend=total_rimr_tend-gr%MS(k,n)%dmassdt(rmt,4)
!!c       end do
!!c
!!c       if(total_mlts_tend>0.0_PS.or.total_rims_tend>0.0_PS) then
!!c          write(*,'("ice>mltck2",3I5,20ES15.6)') KD(n),ID(n),JD(n),&
!!c                   total_mlts_tend*gr%dt,total_mltr_tend*gr%dt,&
!!c                   total_rims_tend*gr%dt,total_rimr_tend*gr%dt
!!c       endif

      icebin_loop2: do i=1,gs%n_bin
        do n=1,gs%L
          total_sloss = &
            ! - collision advection
            max(-gs%MS(i,n)%dmassdt(imt,2),0.0_ds)*real(iupdate_gs(2),PS_KIND) + &
            ! - collision breakup
            ! - riming process with cloud droplets
            max(-gs%MS(i,n)%dmassdt(imt,3),0.0_ds)*real(iupdate_gs(3),PS_KIND) + &
            ! - riming process with rain drops
            max(-gs%MS(i,n)%dmassdt(imt,4),0.0_ds)*real(iupdate_gs(4),PS_KIND) + &
            max(-gs%MS(i,n)%dmassdt(imt,5),0.0_ds)*real(iupdate_gs(5),PS_KIND) + &
            ! - hydrodynamic breakup
            max(-gs%MS(i,n)%dmassdt(imt,6),0.0_ds)*real(iupdate_gs(6),PS_KIND) + &
            ! - secondary nucleation (splintering)
            max(-gs%MS(i,n)%dmassdt(imt,9), 0.0_ds)*real(iupdate_gs(9),PS_KIND) + &
            ! - melting-shedding
            max(-gs%MS(i,n)%dmassdt(imt,10),0.0_ds)*real(iupdate_gs(10),PS_KIND)


          total_sgain = &
            ! - existing mass
            gs%MS(i,n)%mass(imt)/gs%dt + &
!!c               0.9999*gs%MS(i,n)%mass(imt)/gs%dt + &
               ! - collision advection
               max(gs%MS(i,n)%dmassdt(imt,2),0.0_ds)*real(iupdate_gs(2),PS_KIND) + &
!!c               min((rs-gs%MS(i,n)%mass(imt))/gs%dt,max(gs%MS(i,n)%dmassdt(imt,2),0.0_ps)) + &
               ! - riming process with cloud droplets
               max(gs%MS(i,n)%dmassdt(imt,3),0.0_ds)*real(iupdate_gs(3),PS_KIND) + &
               ! - riming process with rain drops
               max(gs%MS(i,n)%dmassdt(imt,4),0.0_ds)*real(iupdate_gs(4),PS_KIND) + &
               ! - collision breakup
               max(gs%MS(i,n)%dmassdt(imt,5),0.0_ds)*real(iupdate_gs(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(gs%MS(i,n)%dmassdt(imt,6),0.0_ds)*real(iupdate_gs(6),PS_KIND) + &
               ! - ice nucleation (contact)
               max(gs%MS(i,n)%dmassdt(imt,8),0.0_ds)*real(iupdate_gs(8),PS_KIND) + &
               ! - secondary nucleation (splintering)
               max( gs%MS(i,n)%dmassdt(imt,9), 0.0_ds)*real(iupdate_gs(9),PS_KIND) + &
               ! - ice nucleation (immersion)
               max( gs%MS(i,n)%dmassdt(imt,11), 0.0_ds)*real(iupdate_gs(11),PS_KIND) + &
               ! - ice nucleation (homogeneous freezing)
               max( gs%MS(i,n)%dmassdt(imt,12), 0.0_ds)*real(iupdate_gs(12),PS_KIND)

!          ierror1(n)=0
          icond1(n)=0
          modc_s(n)=1.0_PS

          if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then

            dum1(1)=max(-gs%MS(i,n)%dmassdt(imt,4),0.0_DS)*real(iupdate_gs(4),PS_KIND)&
                   +max(-gs%MS(i,n)%dmassdt(imt,10),0.0_DS)*real(iupdate_gs(10),PS_KIND)

            modc_s(n) = min(1.0_ps,&
                (total_sgain-dum1(1))/max(total_sloss-dum1(1),total_limit))

            if(total_sloss-dum1(1)<total_limit) then
              modc_s(n)=1.0_PS
            else
              modc_s(n)=max(0.0_PS,modc_s(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_s(n) = mark_acc_s(n) + 1
              icond1(n)=1
            endif
          elseif( total_sgain == 0.0_ps .and. total_sloss == 0.0_ps ) then
             ! this mark should help having bogus tendency in
             ! calculating new tendency (cal_model_tendency)
             gs%MS(i,n)%mark = 4

          endif
          if(modc_s(n)>1.0e+4) then
             LOG_ERROR("cal_mass_budget_col",*) "something wrong",n,gs%MS(i,n)%dmassdt(imt,1:gs%n_tendpros)
             call PRC_abort
          endif
        enddo

        if(all(icond1(1:gs%L)==0)) cycle

        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
          do n=1,gs%L
            icond2(n)=0
            if(gs%MS(i,n)%dmassdt(imt,j) < 0.0_ps) then
              icond2(n)=1
            endif
          enddo
          if(j==2.or.j==5.or.j==6.or.j==9) then
            ! - if those processes contribute to the total loss,
            !   the tendencies in all the bins are reduced for
            !   consistency.
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
             do n = 1, gs%L
             do k = 1, gs%n_bin
              if(icond2(n)==1) then
                gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_s(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_s(n)
              endif
            enddo
            enddo
          else
            ! this may be wrong
            ! write for j=1 same as j=2
            do n=1,gs%L
              if(icond2(n)==1) then
                gs%MS(i,n)%dmassdt(imt,j) = gs%MS(i,n)%dmassdt(imt,j) * modc_s(n)
                acc_mod_s(i,j,n) = acc_mod_s(i,j,n) * modc_s(n)
              end if
            enddo
          endif
        enddo
      enddo icebin_loop2
    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 4. bugdet of each aerosol particles
    if( flagp_a /= 0 ) then
      do ica=1,ncat_a
        aerbin_loop1: do i=1,ga(ica)%n_bin
          do n=1,ga(ica)%L
            total_aloss = &
               ! - ccn activation -
               max(-ga(ica)%MS(i,n)%dmassdt(amt,3),0.0_ds)*real(iupdate_ga(3,ica),PS_KIND)+ &
               ! - ice nucleation (deposition/sorption)
               max(-ga(ica)%MS(i,n)%dmassdt(amt,7),0.0_ds)*real(iupdate_ga(7,ica),PS_KIND)+ &
               ! - ice nucleation (contact)
               max(-ga(ica)%MS(i,n)%dmassdt(amt,8),0.0_ds)*real(iupdate_ga(8,ica),PS_KIND)

!!c             exist_mass=max(ga(ica)%MS(i,n)%mass(1)&
!!c!!c                  -coef4pi3*r3_lmt*ga(ica)%MS(i,n)%den*ga(ica)%MS(i,n)%con,0.0_ps)
!!c                  -coef4pi3*r3_lmt*ga(ica)%MS(i,n)%den*n_lmt_ap,0.0_ps)
            exist_mass=ga(ica)%MS(i,n)%mass(1)

            total_again = &
                  ! - existing mass
!!c                  0.9999*exist_mass/ga(ica)%dt + &
               exist_mass/ga(ica)%dt + &

               ! - vapor deposition (evaporation)
               max(ga(ica)%MS(i,n)%dmassdt(amt,1),0.0_ds)*real(iupdate_ga(1,ica),PS_KIND)

!            ierror1(n)=0
            icond1(n)=0
            modc_a(n)=1.0_PS

            if( total_again < 0.99999_ps*total_aloss.and.total_aloss>total_limit ) then

              dum1(1)=max(-ga(ica)%MS(i,n)%dmassdt(amt,3),0.0_DS)*real(iupdate_ga(3,ica),PS_KIND)&
                  +max(-ga(ica)%MS(i,n)%dmassdt(amt,7),0.0_DS)*real(iupdate_ga(7,ica),PS_KIND)

              modc_a(n) = min(1.0_ps,&
                  (total_again-dum1(1))/max(total_aloss-dum1(1),total_limit))

              if(total_aloss-dum1(1)<total_limit) then
                modc_a(n)=1.0_PS
              else

                modc_a(n)=max(0.0_PS,modc_a(n))

                mark_mod(n) = mark_mod(n) + 1
                mark_acc_a(ica,n) = mark_acc_a(ica,n) + 1

                icond1(n)=1
              endif
            else if( total_again == 0.0_ps .and. total_aloss == 0.0_ps ) then
              ! this mark should help having bogus tendency in
              ! calculating new tendency (cal_model_tendency)
              if( ga(ica)%MS(i,n)%mark == 4 ) then
                ga(ica)%MS(i,n)%mark = 6
              else
                ga(ica)%MS(i,n)%mark = 5
              endif
            endif
            if(modc_a(n)>1.0e+4) then
               LOG_ERROR("cal_mass_budget_col",*) "something wrong",n,ga(ica)%MS(i,n)%dmassdt(amt,1:ga(ica)%n_tendpros)
                call PRC_abort
            endif
            if( (icond1(n) > 0) .and. &
                 ( total_aloss==&
                   max(-ga(ica)%MS(i,n)%dmassdt(amt,3),0.0_DS)*real(iupdate_ga(3,ica),PS_KIND)+&
                   max(-ga(ica)%MS(i,n)%dmassdt(amt,7),0.0_DS)*real(iupdate_ga(7,ica),PS_KIND)) ) then
               LOG_ERROR("cal_mass_budget_col",*) "something is not right at cal_mass_budet a,i,n",i,n
               LOG_ERROR_CONT(*) "total_sloss,sgain",total_aloss,total_again
               do j = 1, ga(ica)%n_tendpros
                  LOG_ERROR("cal_mass_budget_col",*) j,ga(ica)%MS(i,n)%dmassdt(amt,j)
               enddo
               call PRC_abort
            endif
          enddo

          if(all(icond1(1:ga(ica)%L)==0)) cycle

          do j=1,ga(ica)%n_tendpros
            if(iupdate_ga(j,ica)==0) cycle

            do n=1,ga(ica)%L
              icond2(n)=0
              if( ga(ica)%MS(i,n)%dmassdt(amt,j) < 0.0_ps ) then
                icond2(n)=1
              endif
            enddo

            if(j==3) then
              ! ccn activation
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_a(n)
                  gr%ms(k,n)%dmassdt(rmat,j) = gr%ms(k,n)%dmassdt(rmat,j) * modc_a(n)
                  acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gr%L
                mark_acc_r(n) = mark_acc_r(n) + icond1(n)*icond2(n)
              enddo

            elseif(j==7) then
              ! ice nucleation (deposition/sorption)
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_a(n)
                  gs%ms(k,n)%dmassdt(imat,j) = gs%ms(k,n)%dmassdt(imat,j) * modc_a(n)
                  acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gs%L
                mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
              enddo
            elseif(j==8) then
              ! ice nucleation (contact)
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_a(n)
                  acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j)=gr%MS(k,n)%dmassdt(rmt,j) * modc_a(n)
                  acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gr%L
                mark_acc_r(n) = mark_acc_r(n) + icond1(n)*icond2(n)
                mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
              enddo
            endif
          enddo
        enddo aerbin_loop1
      enddo
    endif
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine cal_mass_budget_col

  subroutine cal_mass_budget_vapor( ag, gr, gs, ncat_a,ga, mes_rc, &
       acc_mod_r, acc_mod_s, acc_mod_a,&
       mark_mod, mark_acc_r , mark_acc_s, mark_acc_a, flagp_r, flagp_s, flagp_a,&
       ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    use class_Thermo_Var, only: &
       get_sat_vapor_pres_lk
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate total gain and loss for each group, and modify
    ! only mass tendency
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (AirGroup), intent(inout)    :: ag
    type (Group), intent(inout)         :: gr, gs
    integer,intent(in) :: ncat_a
    type (Group), dimension(ncat_a)          :: ga

    ! message from reality_check subroutine
    ! 0. vapor and hydrometeors do not exist.
    ! 1. vapor only. no hydrometeors.
    ! 2. vapor and liquid hydrometeors exist.
    ! 3. vapor and solid hydrometeoors exist.
    ! 4. vapor, liquid and solid hydrometeors exist.
    !
    integer,dimension(*)   :: mes_rc
    !integer, intent(in)           :: act_type

    integer :: ID(*),JD(*),KD(*)

    real (PS), dimension(mxnbina,mxntend,ncamx,*)    :: acc_mod_a
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_r
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_s
    integer, dimension(*), intent(inout)                :: mark_mod
    integer, dimension(ncamx,*)         :: mark_acc_a
    integer, dimension(*),intent(inout)                :: mark_acc_r, mark_acc_s
    integer, intent(in)                   :: flagp_r, flagp_s, flagp_a
    !
    ! local space
    !
    ! mass of vapor (g/cm^3), saturation mass of vapor
    real(PS),dimension(LMAX)    :: mass_v, mass_vs, mass_vc
    real(PS),dimension(LMAX)    :: rv_new

    real(PS),dimension(LMAX)    :: total_vloss, total_vgain
    real(PS) :: total_rloss, total_rgain, total_sloss, total_sgain,total_aloss, total_again

    ! modification constant for vapor, cloud_drop, rain, and solid_hydro
    real(PS),dimension(LMAX)    :: modc_v, modc_r, modc_s, modc_a
    real(PS),dimension(LMAX)    :: sum_dmassdt_r,sum_dmassdt_s,sum_dmassdt_r_act, &
                                    sum_dmassdt_s_act,sum_dmassdt_s_wet,sum_dmassdt_s_dhf
    integer,dimension(LMAX) :: jice

    integer,dimension(LMAX) :: icond1,icond2,icond3 !,ierror1
    integer,dimension(2,LMAX) :: ick

    ! temperature of freezing
    real(PS), parameter           :: TF = 273.16
    real(PS),parameter :: total_limit=1.0e-30_PS

    real(PS) :: exist_mass

    ! supersaturation at the end of time step
    real(PS),dimension(LMAX) :: svw1,svi1
    real(PS),dimension(2,LMAX) :: rv_sat_n,e_sat_m
    real(PS) :: e_n
    real(PS),dimension(LMAX) :: s_n

    real(PS), parameter :: Rdvchiarui=287.04_PS/461.50_PS

    ! index for deliquescence heterogeneous freezing for ice
    integer,parameter :: jdhf=3
    integer      :: i,j,k, ica,n


    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize
    do n=1,gr%L
      mark_mod(n) = 0


      mass_v(n)=ag%TV(n)%rv*ag%TV(n)%den

      svw1(n)=ag%TV(n)%s_v(1)
      svi1(n)=ag%TV(n)%s_v(2)

      e_sat_m(1,n)=get_sat_vapor_pres_lk(1,ag%TV(n)%T_m,ag%estbar,ag%esitbar)
      e_sat_m(2,n)=get_sat_vapor_pres_lk(2,min(T_0,ag%TV(n)%T_m),ag%estbar,ag%esitbar)
      rv_sat_n(1,n)=Rdvchiarui*e_sat_m(1,n)/(ag%TV(n)%P-e_sat_m(1,n))
      rv_sat_n(2,n)=Rdvchiarui*e_sat_m(2,n)/(ag%TV(n)%P-e_sat_m(2,n))


      sum_dmassdt_r(n)=0.0_PS
      sum_dmassdt_r_act(n)=0.0_PS

      sum_dmassdt_s(n)=0.0_PS
      sum_dmassdt_s_act(n)=0.0_PS
      sum_dmassdt_s_dhf(n)=0.0_PS
      sum_dmassdt_s_wet(n)=0.0_PS

      jice(n)=0
    enddo
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 1. Bugdet of vapor
    !
    !    assumption
    !    - If the obtained positive tendency is too large, then those are
    !      mulplied by some coefficient, which does not vary according to
    !      the size of a hydrometeor in the bin.
    !
    ! NOTE: saturation mixing ratio is recalculated with ambient temperature
    !       because it was calculated for cloud drop surface temperature.

    ! rain group
    do i=1,gr%N_BIN
      do n=1,gr%L
        sum_dmassdt_r(n)=sum_dmassdt_r(n)&
            ! - vapor deposition/evaporation process
            +gr%MS(i,n)%dmassdt(rmt,1)-gr%MS(i,n)%dmassdt(rmat,1)
        sum_dmassdt_r_act(n)=sum_dmassdt_r_act(n)&
            ! - vapor deposition after activation
            +gr%MS(i,n)%dmassdt(rmt,3)-gr%MS(i,n)%dmassdt(rmat,3)
      enddo
    enddo

    ! solid hydrometeor group
    do i=1,gs%N_BIN
      do n=1,gs%L
        sum_dmassdt_s(n)=sum_dmassdt_s(n) &
            ! - vapor deposition/evaporation process
            +gs%MS(i,n)%dmassdt(imt,1)-gs%MS(i,n)%dmassdt(imat,1)
        sum_dmassdt_s_act(n)=sum_dmassdt_s_act(n) &
            ! - vapor deposition after nucleation (deposition nucleation)
            +gs%MS(i,n)%dmassdt(imt,7)-gs%MS(i,n)%dmassdt(imat,7)
        sum_dmassdt_s_dhf(n)=sum_dmassdt_s_dhf(n) &
            ! - deliquence hetero. freezing  process
            +gs%MS(i,n)%dmassdt(imt,3)-gs%MS(i,n)%dmassdt(imat,3)

        ! calculate wet growth of ice
        if(ag%TV(n)%T>=T_0.or.gs%MS(i,n)%inmlt/=0) then
          sum_dmassdt_s_wet(n)=sum_dmassdt_s_wet(n) &
            ! - vapor deposition/evaporation process
!!c            +gs%MS(i,n)%dmassdt(imt,1)-gs%MS(i,n)%dmassdt(imat,1)
            +1.0
        endif
      enddo
    enddo


!    do in=1,gs%n_bin*gs%L
!      n=(in-1)/gs%N_BIN+1
!      i=in-(n-1)*gs%N_BIN
    do n = 1, gs%L
    do i = 1, gs%N_BIN

      if(ag%TV(n)%T<TF) then
        ! --- if ice exists or no hydrometeors ---
        !     the ambient temperature has to be less than freezing, because
        !     otherwise, thin liquid layer produces water saturation condition.
        if(mes_rc(n)>=3) then
          if(gs%MS(i,n)%inmlt==0.and.(gs%MS(i,n)%con>1.0e-30_PS.and.gs%MS(i,n)%mass(imt)>1.0e-30_PS)) then
            jice(n)=1
          endif
        else
        ! - vapor deposition after nucleation (deposition nucleation)
        ! or
        ! - deliquence hetero. freezing
          if(gs%MS(i,n)%dmassdt(imt,7)>1.0e-30 .or.gs%MS(i,n)%dmassdt(imt,3)>1.0e-30) then
            jice(n)=1
          end if
        endif
      endif
    enddo
    enddo

    do n=1,ag%L

      if(jice(n)>0) then
        mass_vs(n)=rv_sat_n(2,n)*ag%TV(n)%den
      else
        mass_vs(n)=rv_sat_n(1,n)*ag%TV(n)%den
      end if

      mass_vc(n)=mass_v(n)

      ! - case of deposition -
      total_vloss(n)=max(sum_dmassdt_r(n),0.0_PS) &
                    +max(sum_dmassdt_s(n),0.0_PS) &
                    +sum_dmassdt_r_act(n)+sum_dmassdt_s_act(n)+sum_dmassdt_s_dhf(n)


      ! - case of evaporation -
      total_vgain(n)=max(-sum_dmassdt_r(n),0.0_PS) &
                    +max(-sum_dmassdt_s(n),0.0_PS)

      total_vgain(n)=total_vgain(n)+mass_vc(n)/gr%dt

      modc_v(n)=1.0_PS
      icond1(n)=0

      if( total_vgain(n) < 0.99999_ps*total_vloss(n).and.total_vloss(n)>total_limit ) then
        modc_v(n) = min(1.0_ps,total_vgain(n)/max(total_vloss(n),total_limit))

        mark_mod(n)=mark_mod(n)+1
        icond1(n)=1

        if(sum_dmassdt_r(n)>0.0_ps) then
          mark_acc_r(n)=mark_acc_r(n)+1
        endif
        if(sum_dmassdt_r_act(n)>0.0_PS) then
          mark_acc_r(n)=mark_acc_r(n)+1
        endif
        if(sum_dmassdt_s(n)>0.0_ps) then
          mark_acc_s(n)=mark_acc_s(n)+1
        endif
        if(sum_dmassdt_s_act(n)>0.0_PS) then
          mark_acc_s(n)=mark_acc_s(n)+1
        endif
        if(sum_dmassdt_s_dhf(n)>0.0_PS) then
          mark_acc_s(n)=mark_acc_s(n)+1
        endif
      endif

    enddo

    ! vapor
    do n=1,ag%L
      if(ag%tv(n)%dmassdt_v(1)>0.0_ps) then
        ag%tv(n)%dmassdt_v(1)=ag%tv(n)%dmassdt_v(1)*modc_v(n)
      end if
      if(ag%tv(n)%dmassdt_v(2)>0.0_ps) then
        ag%tv(n)%dmassdt_v(2)=ag%tv(n)%dmassdt_v(2)*modc_v(n)
      end if
    enddo

    ! aerosol
!    do in=1,ncat_a*ag%L
!      n=(in-1)/ncat_a+1
!      ica=in-(n-1)*ncat_a
    do n = 1, ag%L
    do ica = 1, ncat_a
      if(ica==2.and.&
         sum_dmassdt_s_act(n)>0.0_PS) then

        mark_acc_a(ica,n)=mark_acc_a(ica,n)+icond1(n)
      elseif(sum_dmassdt_r_act(n)>0.0_PS.or.&
             sum_dmassdt_s_dhf(n)>0.0_PS) then
          mark_acc_a(ica,n)=mark_acc_a(ica,n)+icond1(n)
      endif
    enddo
    enddo

    ! rain
!    do in=1,gr%n_bin*gr%L
!      n=(in-1)/gr%N_BIN+1
!      i=in-(n-1)*gr%N_BIN
    do n = 1, gr%L
    do i = 1, gr%N_BIN

      ! - vapor deposition -
      if(sum_dmassdt_r(n)>0.0_PS) then
        gr%MS(i,n)%dmassdt(rmt,1) = gr%MS(i,n)%dmassdt(rmt,1)*modc_v(n)
        gr%MS(i,n)%dmassdt(rmat,1) = gr%MS(i,n)%dmassdt(rmat,1)*modc_v(n)
        acc_mod_r(i,1,n)=acc_mod_r(i,1,n) * modc_v(n)
      endif

      ! - cloud droplets activation process
      if(sum_dmassdt_r_act(n)>0.0_PS) then
        gr%MS(i,n)%dmassdt(rmt,3)= gr%MS(i,n)%dmassdt(rmt,3) * modc_v(n)
        gr%MS(i,n)%dmassdt(rmat,3)= gr%MS(i,n)%dmassdt(rmat,3) * modc_v(n)
        acc_mod_r(i,3,n)=acc_mod_r(i,3,n) * modc_v(n)
      endif
    enddo
    enddo


    ! solid hydrometeor
!    do in=1,gs%n_bin*gs%L
!      n=(in-1)/gs%N_BIN+1
!      i=in-(n-1)*gs%N_BIN
    do n = 1, gs%L
    do i = 1, gs%N_BIN

      ! - vapor deposition -
      if(sum_dmassdt_s(n)>0.0_PS) then

        gs%MS(i,n)%dmassdt(imt,1)=gs%MS(i,n)%dmassdt(imt,1)*modc_v(n)
        gs%MS(i,n)%dmassdt(imat,1)=gs%MS(i,n)%dmassdt(imat,1)*modc_v(n)
        acc_mod_s(i,1,n) = acc_mod_s(i,1,n) * modc_v(n)
      endif

      ! - sorption/deposition ice nucleation process
      if(sum_dmassdt_s_act(n)>0.0_PS) then
        gs%MS(i,n)%dmassdt(imt,7)= gs%MS(i,n)%dmassdt(imt,7) * modc_v(n)
        gs%MS(i,n)%dmassdt(imat,7)= gs%MS(i,n)%dmassdt(imat,7) * modc_v(n)
        acc_mod_s(i,7,n) = acc_mod_s(i,7,n) * modc_v(n)
      endif

      ! - deliquence heterogeneous freezing
      if(sum_dmassdt_s_dhf(n)>0.0_PS) then
        gs%MS(i,n)%dmassdt(imt,3)= gs%MS(i,n)%dmassdt(imt,3) * modc_v(n)
        gs%MS(i,n)%dmassdt(imat,3)= gs%MS(i,n)%dmassdt(imat,3) * modc_v(n)
        acc_mod_s(i,3,n) = acc_mod_s(i,3,n) * modc_v(n)
      endif
    enddo
    enddo

    ! aerosol
    do ica=1,ncat_a
      if(ica==2) then
!        do in=1,ga(ica)%n_bin*ga(ica)%L
!          n=(in-1)/ga(ica)%N_BIN+1
!          i=in-(n-1)*ga(ica)%N_BIN
         do n = 1, ga(ica)%L
         do i = 1, ga(ica)%N_BIN

          if(sum_dmassdt_s_act(n)>0.0_PS) then

            ga(ica)%MS(i,n)%dmassdt(amt,7)= ga(ica)%MS(i,n)%dmassdt(amt,7) * modc_v(n)
            acc_mod_a(i,7,ica,n) = acc_mod_a(i,7,ica,n) * modc_v(n)
          endif

        end do
        end do

      else
!        do in=1,ga(ica)%n_bin*ga(ica)%L
!          n=(in-1)/ga(ica)%N_BIN+1
!          i=in-(n-1)*ga(ica)%N_BIN
         do n = 1, ga(ica)%L
         do i = 1, ga(ica)%N_BIN

          if(sum_dmassdt_r_act(n)>0.0_PS) then

            ! - cloud droplets activation process
            ga(ica)%MS(i,n)%dmassdt(amt,3)= ga(ica)%MS(i,n)%dmassdt(amt,3) * modc_v(n)
            acc_mod_a(i,3,ica,n) = acc_mod_a(i,3,ica,n) * modc_v(n)
          endif

          if(sum_dmassdt_s_dhf(n)>0.0_PS) then

            ! - deliquence hetero freezing -
            ga(ica)%MS(i,n)%dmassdt(amt,6)= ga(ica)%MS(i,n)%dmassdt(amt,6) * modc_v(n)
            acc_mod_a(i,6,ica,n) = acc_mod_a(i,6,ica,n) * modc_v(n)
          endif
        enddo
        enddo
      endif
    enddo


    ! check
    do n=1,ag%L
      if(icond1(n)==1) then
        total_vgain(n)=mass_vc(n)/gr%dt
        total_vloss(n)=0.0_PS

        sum_dmassdt_r(n)=0.0_PS
        sum_dmassdt_r_act(n)=0.0_PS
        sum_dmassdt_s(n)=0.0_PS
        sum_dmassdt_s_act(n)=0.0_PS
        sum_dmassdt_s_dhf(n)=0.0_PS
      endif
    enddo

    ! rain group
    do i=1,gr%N_BIN
      do n=1,gr%L
        if(icond1(n)==1) then
          sum_dmassdt_r(n)=sum_dmassdt_r(n) &
               ! - vapor deposition/evaporation process
               +gr%MS(i,n)%dmassdt(rmt,1)-gr%MS(i,n)%dmassdt(rmat,1)
          sum_dmassdt_r_act(n)=sum_dmassdt_r_act(n) &
               ! - vapor deposition after activation
               +gr%MS(i,n)%dmassdt(rmt,3)-gr%MS(i,n)%dmassdt(rmat,3)
        endif
      end do
    enddo

    ! solid hydrometeor group
    do i=1,gs%N_BIN
      do n=1,gs%L
        if(icond1(n)==1) then
          sum_dmassdt_s(n)=sum_dmassdt_s(n) &
               ! - vapor deposition/evaporation process
               +gs%MS(i,n)%dmassdt(imt,1)-gs%MS(i,n)%dmassdt(imat,1)
          sum_dmassdt_s_act(n)=sum_dmassdt_s_act(n) &
               ! - vapor deposition after nucleation (deposition nucleation)
               +gs%MS(i,n)%dmassdt(imt,7)-gs%MS(i,n)%dmassdt(imat,7)
          sum_dmassdt_s_dhf(n)=sum_dmassdt_s_dhf(n) &
               ! - deliquence hetero. freezing  process
               +gs%MS(i,n)%dmassdt(imt,3)-gs%MS(i,n)%dmassdt(imat,3)
        endif
      enddo
    enddo

    do n=1,ag%L
      if(icond1(n)==1) then
        ! - case of deposition -
        total_vloss(n)=total_vloss(n)+max(sum_dmassdt_r(n),0.0_PS)&
            +max(sum_dmassdt_s(n),0.0_PS)&
            +sum_dmassdt_r_act(n)+sum_dmassdt_s_act(n)+sum_dmassdt_s_dhf(n)

        ! - case of evaporation -
        total_vgain(n)=total_vgain(n)+max(-sum_dmassdt_r(n),0.0_PS)&
            +max(-sum_dmassdt_s(n),0.0_PS)

      endif
    enddo

!    ierror1(1:ag%L)=0
    do n=1,ag%L
      if(mass_v(n)>mass_vs(n) .and.&
         mass_v(n)+total_vgain(n)*gr%dt-min(total_vloss(n),total_vgain(n))*gr%dt<mass_vs(n)) then
         LOG_ERROR("cal_mass_budget_vaporwrite",*) "cal_mass_budget_vapor: too much vapor was used!",KD(n),ID(n),JD(n),&
                     n,mass_v(n),mass_vs(n),total_vgain(n),&
                     total_vloss(n)
          call PRC_abort
      endif
    enddo

    do n=1,ag%L

      ! update the vapor tendency
      ag%tv(n)%dmassdt=max(-mass_v(n)/gr%dt,&
           total_vgain(n)-min(total_vloss(n),total_vgain(n))-mass_vc(n)/gr%dt)

      ! keep the vapor pressure under saturation when vapor produced by evaporation
      rv_new(n)=ag%tv(n)%rv+ag%tv(n)%dmassdt*gr%dt/ag%tv(n)%den

!!c      write(*,*) "ck rv_new",rv_new(n),ag%TV(n)%dmassdt,total_vloss(n),total_vgain(n),&
!!c                  mass_v(n),mass_vc(n)


      e_n=ag%TV(n)%P*rv_new(n)/(Rdvchiarui+rv_new(n))
      s_n(n)=e_n/e_sat_m(1,n)-1.0

      modc_v(n)=1.0_PS
      icond1(n)=0
!      ierror1(n)=0

      if(mes_rc(n)==3.and.abs(sum_dmassdt_s_wet(n))<=1.0e-25_PS) then
        ! this equilibrium condition over ice saturation in case of evaporation
        ! should occur when only ice particles exist. If liquid hydrometeors
        ! coexists, and they evaporate, then still supersaturation over ice
        ! should be possible.

        if((rv_new(n)-rv_sat_n(2,n))/rv_sat_n(2,n)>1.0e-4.and.&
           (ag%tv(n)%rv<ag%tv(n)%rv_sat(2).or.svi1(n)<0.0_ps)) then

          if(max(0.0_ps,-sum_dmassdt_r(n))+max(0.0_ps,-sum_dmassdt_s(n))>total_limit) then

            icond1(n)=1
            modc_v(n)=max(0.0_PS,(rv_sat_n(2,n)-ag%tv(n)%rv)*ag%tv(n)%den/gr%dt/&
                  (max(0.0_ps,-sum_dmassdt_r(n))+max(0.0_ps,-sum_dmassdt_s(n))))

          endif

        endif
      else
        if( (rv_new(n)-rv_sat_n(1,n))/rv_sat_n(1,n)>1.0e-4_PS.and.&
             ((ag%tv(n)%rv<ag%tv(n)%rv_sat(1).or.svw1(n)<0.0_ps).or.&
               (ag%tv(n)%s_v(1)>0.0_PS.and.s_n(n)<0.0_PS).or.&
              gr%mark_er(n)==1)) then

          if(max(0.0_ps,-sum_dmassdt_r(n))+max(0.0_ps,-sum_dmassdt_s(n))>total_limit) then
!!c             modc_v=(ag%tv(n)%rv_sat(1)-ag%tv(n)%rv)*ag%tv(n)%den/gr%dt/&

             icond1(n)=2
             modc_v(n)=max(0.0_PS,(rv_sat_n(1,n)-ag%tv(n)%rv)*ag%tv(n)%den/gr%dt/&
                         (max(0.0_ps,-sum_dmassdt_r(n))+max(0.0_ps,-sum_dmassdt_s(n))))
          endif
        endif
      endif

      if(modc_v(n)>=0.0_ps.and.modc_v(n)<0.99999_PS) then
        if(sum_dmassdt_r(n)<0.0_ps) then
          ! liquid hydrometeors
          mark_acc_r(n)=mark_acc_r(n)+1
        endif
        if(sum_dmassdt_s(n)<0.0_ps) then
          ! solid hydrometeor group
          mark_acc_s(n) = mark_acc_s(n) + 1
        endif
      elseif(modc_v(n)>1.0) then
         LOG_ERROR("cal_mass_budget_vapor",*) "k,i,j,n",KD(n),ID(n),JD(n),n
         LOG_ERROR_CONT(*) "sum_dmassdt_r,vc",sum_dmassdt_r(n),mass_vc(n)
         LOG_ERROR_CONT(*) "dmassdt_v",ag%tv(n)%dmassdt
         LOG_ERROR_CONT(*) "gain,loss",total_vgain(n),total_vloss(n)

          total_vgain(n)=mass_vc(n)/gr%dt
          total_vloss(n)=0.0_ps
          ! rain group
          sum_dmassdt_r(n)=0.0
          do i = 1, gr%n_bin
             ! - vapor deposition/evaporation process
             sum_dmassdt_r(n)=sum_dmassdt_r(n)+gr%MS(i,n)%dmassdt(rmt,1)&
               +max(0.0_ds,gr%MS(i,n)%dmassdt(rmt,3)-gr%MS(i,n)%dmassdt(rmat,3))
          end do
          ! - case of deposition -
          total_vloss(n)=total_vloss(n)+max(sum_dmassdt_r(n),0.0_ps)&
               +max(sum_dmassdt_s(n),0.0_ps)

          ! - case of evaporation -
          total_vgain(n)=total_vgain(n)+max(-sum_dmassdt_r(n),0.0_ps)&
               +max(-sum_dmassdt_s(n),0.0_ps)

          LOG_ERROR_CONT(*) "loss,gain,mod",total_vloss(n),total_vgain(n),modc_v(n)
          LOG_ERROR_CONT('("dmassdt 1",25es15.6)') (gr%MS(i,n)%dmassdt(rmt,1),i=1,gr%n_bin)
          LOG_ERROR_CONT('("dmassdt 3",25es15.6)') (gr%MS(i,n)%dmassdt(rmt,3),i=1,gr%n_bin)
          call PRC_abort
      endif
    enddo

    do n=1,gr%L
      icond2(n)=0
      icond3(n)=0
      if(sum_dmassdt_r(n)<0.0_ps) then
        icond2(n)=1
      endif
      if(sum_dmassdt_s(n)<0.0_ps) then
        icond3(n)=1
      endif
    enddo

    do n=1,gr%L
      if(icond2(n)==1) then
        ag%tv(n)%dmassdt_v(1)=ag%tv(n)%dmassdt_v(1)*modc_v(n)
      endif
      if(icond3(n)==1) then
        ag%tv(n)%dmassdt_v(2)=ag%tv(n)%dmassdt_v(2)*modc_v(n)
      endif
    enddo

    ! - vapor deposition/evaporation process
!    do in=1,gr%n_bin*gr%L
!      n=(in-1)/gr%N_BIN+1
!      i=in-(n-1)*gr%N_BIN
    do n = 1, gr%L
    do i = 1, gr%N_BIN

      if(icond2(n)==1) then
        gr%MS(i,n)%dmassdt(rmt,1)=gr%MS(i,n)%dmassdt(rmt,1)*modc_v(n)
        gr%MS(i,n)%dmassdt(rmat,1)=gr%MS(i,n)%dmassdt(rmat,1)*modc_v(n)
        acc_mod_r(i,1,n) = acc_mod_r(i,1,n) * modc_v(n)

      endif
    enddo
    enddo

    do n=1,gr%L
      ick(1,n)=0
      ick(2,n)=0
    enddo
!    do in=1,gr%n_bin*gr%L
!      n=(in-1)/gr%N_BIN+1
!      i=in-(n-1)*gr%N_BIN
    do n = 1, gr%L
    do i = 1, gr%N_BIN
      if(icond2(n)==1.and.gr%MS(i,n)%inevp>0) then
        ick(gr%MS(i,n)%inevp,n)=1
      endif
    enddo
    enddo

!    do in=1,2*gr%L
!      n=(in-1)/2+1
!      ica=in-(n-1)*2
    do n = 1, gr%L
    do ica = 1, 2
      if(ick(ica,n)>0) then
        ga(ica)%ms(1,n)%dmassdt(amt,1)= ga(ica)%ms(1,n)%dmassdt(amt,1) * modc_v(n)
        acc_mod_a(1,1,ica,n) = acc_mod_a(1,1,ica,n) * modc_v(n)

        mark_acc_a(ica,n)=mark_acc_a(ica,n) + min(icond1(n),1)
      endif
    enddo
    enddo

    ! - vapor deposition/evaporation process
!    do in=1,gs%n_bin*gs%L
!      n=(in-1)/gs%N_BIN+1
!      i=in-(n-1)*gs%N_BIN
    do n = 1, gs%L
    do i = 1, gs%N_BIN
      if(icond3(n)==1) then
        gs%MS(i,n)%dmassdt(imt,1)=gs%MS(i,n)%dmassdt(imt,1)*modc_v(n)
        gs%MS(i,n)%dmassdt(imat,1)=gs%MS(i,n)%dmassdt(imat,1)*modc_v(n)
        acc_mod_s(i,1,n) = acc_mod_s(i,1,n) * modc_v(n)
      endif
    enddo
    enddo

    do n=1,gr%L
      ick(1,n)=0
      ick(2,n)=0
    enddo
!    do in=1,gs%n_bin*gs%L
!      n=(in-1)/gs%N_BIN+1
!      i=in-(n-1)*gs%N_BIN
    do n = 1, gs%L
    do i = 1, gs%N_BIN
      if(icond3(n)==1.and.gs%MS(i,n)%inevp>0) then
        ick(gs%MS(i,n)%inevp,n)=1
      endif
    enddo
    enddo
    do ica=1,2
      do n=1,gs%L
        if(ick(ica,n)>0) then
          ga(ica)%ms(1,n)%dmassdt(amt,1)= ga(ica)%ms(1,n)%dmassdt(amt,1) * modc_v(n)
          acc_mod_a(1,1,ica,n) = acc_mod_a(1,1,ica,n) * modc_v(n)

          mark_acc_a(ica,n)=mark_acc_a(ica,n) + min(icond1(n),1)
        endif
      enddo
    enddo

    do n=1,ag%L
!      ierror1(n)=0

      if(icond1(n)==1) then
        ag%tv(n)%dmassdt=(rv_sat_n(2,n)-ag%tv(n)%rv)/gr%dt*ag%tv(n)%den
      elseif(icond1(n)==2) then
        ag%tv(n)%dmassdt=(rv_sat_n(1,n)-ag%tv(n)%rv)/gr%dt*ag%tv(n)%den
      endif
      if(ag%tv(n)%dmassdt>1.0e+2_PS.or.ag%tv(n)%dmassdt<-1.0e+2_PS) then
         LOG_ERROR("cal_mass_budget_vapor",*) "cal_mass_budget > vp dmdt is out of reality"
         LOG_ERROR_CONT(*) "dmdt,tot_g,tot_l",n,ag%tv(n)%dmassdt,total_vgain,total_vloss
         LOG_ERROR_CONT(*) "sum",sum_dmassdt_r(n),sum_dmassdt_s(n),sum_dmassdt_r_act(n), &
              sum_dmassdt_s_act(n),sum_dmassdt_s_dhf(n)
         call PRC_abort
      endif

      rv_new(n)=ag%tv(n)%rv+ag%tv(n)%dmassdt*gr%dt/ag%tv(n)%den
      e_n=ag%TV(n)%P*rv_new(n)/(Rdvchiarui+rv_new(n))
      s_n(n)=e_n/e_sat_m(1,n)-1.0_PS
    enddo

    if ( debug ) then
       do n=1,ag%L
!      ierror1(n)=0
 !kid    if(s_n>0.1) then
          if(s_n(n)>0.5_PS) then
            write(*,*) "cal_mass_budget_vapor > super saturated at",KD(n),ID(n),JD(n),jice(n)
            write(*,*) "n rv,rvs",rv_new(n),ag%tv(n)%rv_sat(1),rv_sat_n(1,n),sum_dmassdt_s_wet(n)
            write(*,*) "T, sw,si,sw2,si2,s_n,T_m,w",ag%TV(n)%T,ag%tv(n)%s_v(1),ag%tv(n)%s_v(2),&
                 ag%tv(n)%s_v_n(1),ag%tv(n)%s_v_n(2),s_n(n),ag%TV(n)%T_m,ag%tv(n)%w

            total_vgain(n)=mass_vc(n)/gr%dt
            total_vloss(n)=0.0_PS
            sum_dmassdt_r(n)=0.0_PS
            sum_dmassdt_r_act(n)=0.0_PS
            sum_dmassdt_s(n)=0.0_PS
            sum_dmassdt_s_act(n)=0.0_PS
            sum_dmassdt_s_dhf(n)=0.0_PS

            ! rain group
            sum_dmassdt_r(n)=sum_dmassdt_r(n) &
                 ! - vapor deposition/evaporation process
                 +gr%MS(i,n)%dmassdt(rmt,1)-gr%MS(i,n)%dmassdt(rmat,1)
            sum_dmassdt_r_act(n)=sum_dmassdt_r_act(n) &
                 ! - vapor deposition after activation
                 +max(0.0_DS,gr%MS(i,n)%dmassdt(rmt,3)-gr%MS(i,n)%dmassdt(rmat,3))

            ! solid hydrometeor group
            sum_dmassdt_s(n)=sum_dmassdt_s(n) &
                 ! - vapor deposition/evaporation process
                 +gs%MS(i,n)%dmassdt(imt,1)-gs%MS(i,n)%dmassdt(imat,1)
            sum_dmassdt_s_act(n)=sum_dmassdt_s_act(n) &
                 ! - vapor deposition after nucleation (deposition nucleation)
                 +max(0.0_DS,gs%MS(i,n)%dmassdt(imt,7)-gs%MS(i,n)%dmassdt(imat,7))
            sum_dmassdt_s_dhf(n)=sum_dmassdt_s_dhf(n) &
                 ! - deliquence hetero. freezing  process
                 +max(0.0_DS,gs%MS(i,n)%dmassdt(imt,3)-gs%MS(i,n)%dmassdt(imat,3))

            ! - case of deposition -
            total_vloss(n)=total_vloss(n)+max(sum_dmassdt_r(n),0.0_PS)&
                 +max(sum_dmassdt_s(n),0.0_PS)&
                 +sum_dmassdt_r_act(n)+sum_dmassdt_s_act(n)

            ! - case of evaporation -
            total_vgain(n)=total_vgain(n)+max(-sum_dmassdt_r(n),0.0_PS)&
                 +max(-sum_dmassdt_s(n),0.0_PS)

            write(*,*) "n",n
            write(*,*) "loss,gain,mod",total_vloss(n),total_vgain(n),modc_v(n)
            write(*,*) "dmassdt 1", (gr%MS(i,n)%dmassdt(rmt,1),i=1,gr%n_bin)
            write(*,*) "dmassdt 3", (gr%MS(i,n)%dmassdt(rmt,3),i=1,gr%n_bin)
            write(*,*) "rain", (gr%MS(i,n)%mass(rmt),i=1,gr%n_bin)
            write(*,*) "dmassdt 1 sol", (gs%MS(i,n)%dmassdt(imt,1),i=1,gs%n_bin)
            write(*,*) "solid", (gs%MS(i,n)%mass(imt),i=1,gs%n_bin)
            write(*,*) "solid tmp", (gs%MS(i,n)%tmp,i=1,gs%n_bin)
            write(*,*) "solid inmlt", (gs%MS(i,n)%inmlt,i=1,gs%n_bin)
            write(*,*) "ap con,mass 1", ga(1)%ms(1,n)%con,ga(1)%ms(1,n)%mass(1)

            write(*,*) "massvc,depr,actr,depi,acti,dhfi",mass_vc(n)/gr%dt,&
                 sum_dmassdt_r(n),sum_dmassdt_r_act(n),&
                 sum_dmassdt_s(n),sum_dmassdt_s_act(n),sum_dmassdt_s_dhf(n)
            write(*,*) " "
          end if
       enddo
    endif


    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2. budget of vapor-related processes
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2.1. bugdet of each rain bin
    ! calculate total tendency of riming of rain drops by ices
    if( flagp_r > 0 ) then

      liqbin_loop1: do i=1,gr%n_bin

        do n=1,gr%L

          total_rloss = &
               ! - evaporation
               max(-gr%MS(i,n)%dmassdt(rmt,1),0.0_ds)

          total_rgain = &
               ! - existing mass
               gr%MS(i,n)%mass(rmt)/gr%dt + &
               ! - vapor deposition
               max(gr%MS(i,n)%dmassdt(rmt,1),0.0_ds) + &
               ! - production of cloud droplets by ccn activation
               max(gr%MS(i,n)%dmassdt(rmt,3),0.0_ds)

!          ierror1(n)=0
          icond1(n)=0
          modc_r(n)=1.0_PS

          if( total_rgain < 0.99999_ps*total_rloss.and.total_rloss>total_limit) then

            modc_r(n) = min(1.0_ps,total_rgain/max(total_rloss,total_limit))

            mark_mod(n) = mark_mod(n) + 1
            mark_acc_r(n) = mark_acc_r(n) + 1
            icond1(n)=1

!!c             write(*,22) KD(n),ID(n),JD(n),i,modc_r,total_rgain,total_rloss
22           format("mass_vapor:modc_r,rgain,rloss",4i5,3es15.6)
          elseif( total_rgain == 0.0_ps .and. total_rloss == 0.0_ps ) then
            ! this mark should help having bogus tendency in
            ! calculating new tendency (cal_model_tendency)
            gr%MS(i,n)%mark = 4
          endif
          if(modc_r(n)>1.0e+4) then
             LOG_ERROR("cal_mass_budget_vapor",*) "something wrong",n,gr%MS(i,n)%dmassdt(rmt,1:gr%n_tendpros)
             call PRC_abort
         endif
        enddo

        if(all(icond1(1:gr%L)==0)) cycle

        do j=1,gr%n_tendpros
          if(j/=1.and.j/=3) cycle

          do n=1,gr%L
            icond2(n)=0
            if( gr%MS(i,n)%dmassdt(rmt,j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo

          if(j==1) then
            ! evaporation
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if(icond2(n)==1) then
                gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_r(n)
                gr%MS(k,n)%dmassdt(rmat,j) = gr%MS(k,n)%dmassdt(rmat,j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n=1,gr%L
              ick(1,n)=0
              if(icond2(n)==1.and.gr%MS(i,n)%inevp>0) then
                ick(1,n)=gr%MS(i,n)%inevp
              endif
            enddo

            do n=1,gr%L
              if(ick(1,n)>0) then
                ica=ick(1,n)

                ga(ica)%ms(1,n)%dmassdt(amt,1)= ga(ica)%ms(1,n)%dmassdt(amt,1) * modc_r(n)
                acc_mod_a(1,1,ica,n) = acc_mod_a(1,1,ica,n) * modc_r(n)

                mark_acc_a(ica,n)=mark_acc_a(ica,n) + icond1(n)
              endif
            enddo
          else
            do n=1,gr%L
              if(icond2(n)==1) then
                gr%MS(i,n)%dmassdt(rmt,j) = gr%MS(i,n)%dmassdt(rmt,j) * modc_r(n)
                acc_mod_r(i,j,n) = acc_mod_r(i,j,n) * modc_r(n)
              endif
            enddo
          endif
        enddo
      enddo liqbin_loop1
    end if
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2.2. bugdet of each solid hydrometeor bin
    if( flagp_s > 0 ) then
      icebin_loop1: do i=1,gs%n_bin

        do n=1,gs%L

          total_sloss = &
               ! - evaporation
               max(-gs%MS(i,n)%dmassdt(imt,1),0.0_ds)


          total_sgain = &
               ! - existing mass
               gs%MS(i,n)%mass(imt)/gs%dt + &
               ! - vapor deposition
               max(gs%MS(i,n)%dmassdt(imt,1),0.0_ds) + &
               ! - ice nucleation (sorption/deposition)
               max(gs%MS(i,n)%dmassdt(imt,7),0.0_ds) + &
               ! - ice nucleation (deliquence-heterogeneous freezing)
               max(gs%MS(i,n)%dmassdt(imt,3),0.0_ds)

!          ierror1(n)=0
          icond1(n)=0
          modc_s(n)=1.0_PS

          if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then

            modc_s(n) = min(1.0_ps,total_sgain/max(total_sloss,total_limit))

            mark_mod(n) = mark_mod(n) + 1
            mark_acc_s(n) = mark_acc_s(n) + 1

            icond1(n)=1

!!c             write(*,23) i,n,modc_s,total_sgain,total_sloss
23           format("mass_vapor:modc_s,sgain,sloss",2i5,3es15.6)
          else if( total_sgain == 0.0_ps .and. total_sloss == 0.0_ps ) then
             ! this mark should help having bogus tendency in
             ! calculating new tendency (cal_model_tendency)
             gs%MS(i,n)%mark = 4

          end if
          if(modc_s(n)>1.0e+4) then
             LOG_ERROR("cal_mass_budget_vapor",*) "something wrong",n,gs%MS(i,n)%dmassdt(imt,1:gs%n_tendpros)
             call PRC_abort
          endif
        enddo

        if(all(icond1(1:gs%L)==0)) cycle

        do j=1,gs%n_tendpros
          if(j/=1.and.j/=7.and.j/=3) cycle

          do n=1,gs%L
            icond2(n)=0
            if( gs%MS(i,n)%dmassdt(imt,j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo

          if(j==1) then
            ! evaporation
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
             do n = 1, gs%L
             do k = 1, gs%n_bin

              if( icond2(n)==1) then
                gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_s(n)
                gs%MS(k,n)%dmassdt(imat,j) = gs%MS(k,n)%dmassdt(imat,j) * modc_s(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_s(n)
              endif
            enddo
            enddo
            do n=1,gs%L
              ick(1,n)=0
              if(icond2(n)==1.and.gs%MS(i,n)%inevp>0) then
                ick(1,n)=gs%MS(i,n)%inevp
              endif
            enddo
            do n=1,gs%L
              if(ick(1,n)>0) then
                ica=ick(1,n)
                ga(ica)%ms(1,n)%dmassdt(amt,1)= ga(ica)%ms(1,n)%dmassdt(amt,1) * modc_s(n)
                acc_mod_a(1,1,ica,n) = acc_mod_a(1,1,ica,n) * modc_s(n)

                mark_acc_a(ica,n)=mark_acc_a(ica,n) + icond1(n)
              endif
            enddo
          else
            ! this may be wrong
            ! write for j=1 same as j=2
            do n=1,gs%L
              if( icond2(n)==1) then
                gs%MS(i,n)%dmassdt(imt,j) = gs%MS(i,n)%dmassdt(imt,j) * modc_s(n)
                acc_mod_s(i,j,n) = acc_mod_s(i,j,n) * modc_s(n)
              endif
            end do
          endif
        enddo
      enddo icebin_loop1
    endif

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 4. bugdet of each aerosol particles
    if( flagp_a /= 0 ) then
      do ica=1,ncat_a
        aerbin_loop1: do i=1, ga(ica)%n_bin
          do n=1,ag%L

            total_aloss = &
                  ! - ccn activation -
                  max(-ga(ica)%MS(i,n)%dmassdt(amt,3),0.0_ds)+ &
                  ! - ice nucleation (deposition/sorption)
                  max(-ga(ica)%MS(i,n)%dmassdt(amt,7),0.0_ds)+ &
                  ! - ice nucleation (deliquence heterogeneous freezing)
                  max(-ga(ica)%MS(i,n)%dmassdt(amt,6),0.0_ds)

            exist_mass=ga(ica)%MS(i,n)%mass(1)

            total_again = &
                  ! - existing mass
                  exist_mass/ga(ica)%dt + &
                  ! - vapor deposition (evaporation)
                  max(ga(ica)%MS(i,n)%dmassdt(amt,1),0.0_ds)

!            ierror1(n)=0
            icond1(n)=0
            modc_a(n)=1.0_PS

            if( total_again < 0.99999_ps*total_aloss.and.total_aloss>total_limit ) then
              modc_a(n) = min(1.0_ps,total_again/max(total_aloss,total_limit))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_a(ica,n) = mark_acc_a(ica,n) + 1

              icond1(n)=1

            else if( total_again == 0.0_ps .and. total_aloss == 0.0_ps ) then
              ! this mark should help having bogus tendency in
              ! calculating new tendency (cal_model_tendency)
              if( ga(ica)%MS(i,n)%mark == 4 ) then
                ga(ica)%MS(i,n)%mark = 6
              else
                ga(ica)%MS(i,n)%mark = 5
              end if
            endif
            if(modc_a(n)>1.0e+4) then
               LOG_ERROR("cal_mass_budget_vapor",*) "something wrong",ga(ica)%MS(i,n)%dmassdt(amt,1:ga(ica)%n_tendpros)
                call PRC_abort
            end if
          enddo

          if(all(icond1(1:ga(ica)%L)==0)) cycle

          do j=1,ga(ica)%n_tendpros
            if(j/=1.and.j/=3.and.j/=7.and.j/=6) cycle

            do n=1,gr%L
              icond2(n)=0
              if( ga(ica)%MS(i,n)%dmassdt(amt,j) < 0.0_ps ) then
                icond2(n)=1
              endif
            enddo

            if(j==3) then
              ! ccn activation
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if( icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j) = gr%MS(k,n)%dmassdt(rmt,j) * modc_a(n)
                  gr%ms(k,n)%dmassdt(rmat,j) = gr%ms(k,n)%dmassdt(rmat,j) * modc_a(n)
                  acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_a(n)
                endif
              end do
              end do
              do n=1,gr%L
                mark_acc_r(n)=mark_acc_r(n)+icond1(n)*icond2(n)
              enddo

            elseif(j==7) then
              ! ice nucleation (deposition/sorption)
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_a(n)
                  gs%ms(k,n)%dmassdt(imat,j) = gs%ms(k,n)%dmassdt(imat,j) * modc_a(n)
                  acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gs%L
                mark_acc_s(n)=mark_acc_s(n)+icond1(n)*icond2(n)
              enddo
            elseif(j==6) then
              ! ice nucleation (deliquence heterogeneous freezing)
              !
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%MS(k,n)%dmassdt(imt,jdhf) = gs%MS(k,n)%dmassdt(imt,jdhf) * modc_a(n)
                  gs%ms(k,n)%dmassdt(imat,jdhf) = gs%ms(k,n)%dmassdt(imat,jdhf) * modc_a(n)
                  acc_mod_s(k,jdhf,n) = acc_mod_s(k,jdhf,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gs%L
                mark_acc_s(n)=mark_acc_s(n)+icond1(n)*icond2(n)
              enddo
            elseif(j==8) then
              ! ice nucleation (contact)
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dmassdt(amt,j) = ga(ica)%ms(k,n)%dmassdt(amt,j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%MS(k,n)%dmassdt(imt,j) = gs%MS(k,n)%dmassdt(imt,j) * modc_a(n)
                  gs%ms(k,n)%dmassdt(imat,j) = gs%ms(k,n)%dmassdt(imat,j) * modc_a(n)
                  acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j)=gr%MS(k,n)%dmassdt(rmt,j) * modc_a(n)
                  acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,ag%L
                mark_acc_s(n)=mark_acc_s(n)+icond1(n)*icond2(n)
                mark_acc_r(n)=mark_acc_r(n)+icond1(n)*icond2(n)
              enddo
            endif
          enddo
        enddo aerbin_loop1
      enddo
    endif
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine cal_mass_budget_vapor


  subroutine cal_con_budget_col(gr, gs, ncat_a,ga, &
       acc_mod_r, acc_mod_s,acc_mod_a, &
       mark_mod, mark_acc_r , mark_acc_s, mark_acc_a, flagp_r, flagp_s, flagp_a,&
       iupdate_gr,iupdate_gs,iupdate_ga&
       )
    use scale_prc, only: &
       PRC_abort
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate total gain and loss for each group, and modify
    ! only concentration tendency
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (group), intent(inout)         :: gr, gs
    integer,intent(in) :: ncat_a
    type (Group), dimension(ncat_a)          :: ga
    !
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_r
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_s
    real (PS), dimension(mxnbina,mxntend,ncamx,*)    :: acc_mod_a
    !
    integer, dimension(*), intent(inout)     :: mark_mod
    integer, dimension(ncamx,*)    :: mark_acc_a
    integer, dimension(*),intent(inout)     :: mark_acc_r, mark_acc_s
    integer, intent(in)                   :: flagp_a, flagp_r, flagp_s
    !integer :: ID(*),JD(*),KD(*)
    integer,dimension(*),intent(in) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,*),intent(in) :: iupdate_ga
    !
    ! local space
    !
    real(ps),dimension(LMAX)    :: total_rimr_tend, total_rims_tend,total_mlts_tend,total_mltr_tend
    real (ps)    :: total_aloss, total_again, &
         total_rloss, total_rgain, total_sloss, total_sgain
    ! modification constant for vapor, cloud_drop, rain, and solid_hydro
    real(ps),dimension(LMAX)    :: modc_a, modc_r, modc_s
    real(ps)    :: dum1
    real(ps),dimension(LMAX)    :: mod_dum
    integer,dimension(LMAX) :: ick
    integer,dimension(LMAX) :: icond1,icond2 !,ierror1
    real(ps),parameter :: total_limit=1.0e-30_ps
    integer      :: i,j,k,ica,n

    !     the minimum possible background concentration for large particle
    !     is assumed to be 1.0e-5 cm^-3
    real(ps),parameter :: n_lmt_ap=1.0e-5
!parcel model    real(ps),parameter :: n_lmt_ap=1.0e-15
    real(ps) :: exist_con

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize
    do n=1,gr%L
      mark_mod(n) = 0
    enddo
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2. bugdet of each rain bin
    ! calculate total tendency of riming of rain drops by ices
    if( flagp_r > 0 ) then

      liqbin_loop1: do i=1,gr%n_bin
        do n=1,gr%L

          total_rloss = &
               ! - collision advection
               max(-gr%MS(i,n)%dcondt(2),0.0_ds)*real(iupdate_gr(2),PS_KIND) + &
               ! - riming process by a solid hyrometeor
               max(-gr%MS(i,n)%dcondt(4), 0.0_ds)*real(iupdate_gr(4),PS_KIND) + &
               ! - collision breakup
               max(-gr%MS(i,n)%dcondt(5),0.0_ds)*real(iupdate_gr(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(-gr%MS(i,n)%dcondt(6),0.0_ds)*real(iupdate_gr(6),PS_KIND) + &
               ! - loss of cloud droplets by contact nucleation
               max(-gr%MS(i,n)%dcondt(8),0.0_ds)*real(iupdate_gr(8),PS_KIND) + &
               ! - auto conversion
               max(-gr%MS(i,n)%dcondt(9), 0.0_ds)*real(iupdate_gr(9),PS_KIND) + &
               ! - ice nucleation (immersion)
               max(-gr%MS(i,n)%dcondt(11), 0.0_ds)*real(iupdate_gr(11),PS_KIND) + &
               ! - ice nucleation (homogeneous freezing)
               max(-gr%MS(i,n)%dcondt(12), 0.0_ds)*real(iupdate_gr(12),PS_KIND)

          total_rgain = &
               ! - existing mass
               gr%MS(i,n)%con/gr%dt + &
!!c               0.9999*gr%MS(i,n)%con/gr%dt + &
               ! - collision advection
               max(gr%MS(i,n)%dcondt(2),0.0_ds)*real(iupdate_gr(2),PS_KIND) + &
               ! - collision breakup
               max(gr%MS(i,n)%dcondt(5),0.0_ds)*real(iupdate_gr(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(gr%MS(i,n)%dcondt(6),0.0_ds)*real(iupdate_gr(6),PS_KIND) + &
!!c               ! - auto-conversion
!!c               max(gr%MS(i,n)%dcondt(7),0.0_ds) + &
               ! - auto conversion
               max(gr%MS(i,n)%dcondt(9), 0.0_ds)*real(iupdate_gr(9),PS_KIND) + &
               ! - melting-shedding
               max(gr%MS(i,n)%dcondt(10),0.0_ds)*real(iupdate_gr(10),PS_KIND)

!          ierror1(n)=0
          icond1(n)=0
          modc_r(n)=1.0_PS

          if( total_rgain < 0.99999_ps*total_rloss.and.total_rloss>total_limit) then

            if( total_rloss==max(-gr%MS(i,n)%dcondt(1),0.0_DS)) then
!!c                write(*,*) "something is not right at cal_con_budet"
!!c                write(*,*) "total_rloss,rgain",total_rloss,total_rgain
!!c                do j = 1, gr%n_tendpros
!!c                   write(*,*) j,gr%MS(i,n)%dcondt(j)
!!c                end do
!!c                stop
              gr%MS(i,n)%dcondt(1)=-total_rgain
            else

              modc_r(n) = min(1.0_ps,&
                      total_rgain/max(total_rloss,total_limit))

              modc_r(n)=max(0.0_PS,modc_r(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_r(n) = mark_acc_r(n) + 1

              icond1(n)=1
            end if
          else if( total_rgain == 0.0_ps .and. total_rloss == 0.0_ps ) then
            ! this mark should help having bogus tendency in
            ! calculating new tendency (cal_model_tendency)
            if( gr%MS(i,n)%mark == 4 ) then
              gr%MS(i,n)%mark = 6
            else
              gr%MS(i,n)%mark = 5
            end if
          endif
          if(modc_r(n)>1.0e+4) then
             LOG_ERROR("cal_con_budget_col",*) "something wrong",n,gr%MS(i,n)%dcondt(1:gr%n_tendpros)
             call PRC_abort
          endif
        enddo

        if(all(icond1(1:gr%L)==0)) cycle

        do j=1,gr%n_tendpros
          if(iupdate_gr(j)==0) cycle

          do n=1,gr%L
            icond2(n)=0
            if( gr%MS(i,n)%dcondt(j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo

          if(j == 2 .or. j == 5 .or. j == 6.or.j==9) then
            ! - if those processes contribute to the total loss,
            !   the tendencies in all the bins are reduced for
            !   consistency.
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if(icond2(n)==1) then
                gr%ms(k,n)%dcondt(j) = gr%ms(k,n)%dcondt(j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            end do
            end do
          elseif(j==4.and.flagp_s>0) then
            ! riming process
            do n=1,gr%L
              if( icond2(n)==1) then
                gr%MS(i,n)%dcondt(j)=gr%MS(i,n)%dcondt(j)*modc_r(n)
                acc_mod_r(i,j,n)=acc_mod_r(i,j,n)*modc_r(n)
              endif
              total_rimr_tend(n)=0.0_ps
              total_rims_tend(n)=0.0_ps
            enddo
            do k=1,gr%n_bin
              do n=1,gr%L
                total_rimr_tend(n)=total_rimr_tend(n)-gr%MS(k,n)%dmassdt(rmt,j)*acc_mod_r(k,j,n)
              enddo
            enddo
            do k=1,gs%n_bin
              do n=1,gr%L
                total_rims_tend(n)=total_rims_tend(n)+gs%MS(k,n)%dmassdt(imt,j)*acc_mod_s(k,j,n)
              enddo
            enddo
            do n=1,gr%L
              mod_dum(n)=1.0_PS
              if( icond2(n)==1) then
                mod_dum(n)=max(0.0_PS,min(1.0_PS,total_rimr_tend(n)/max(total_limit,total_rims_tend(n))))
              endif
              mark_acc_s(n)=mark_acc_s(n)+icond1(n)*icond2(n)
            enddo
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
            do n = 1, gs%L
            do k = 1, gs%n_bin
              gs%ms(k,n)%dcondt(j)=gs%ms(k,n)%dcondt(j)*mod_dum(n)
              acc_mod_s(k,j,n)=acc_mod_s(k,j,n)*mod_dum(n)
            enddo
            enddo

!!c                      write(*,'("rain con>rim",4I5,10ES15.6)') KD(n),ID(n),JD(n),i,&
!!c                             mod_dum,modc_r,total_rims_tend*gr%dt,total_rimr_tend*gr%dt

          else if((j == 11 .or. j == 12).and. flagp_s > 0 ) then
            ! - immersion freezing or homogeneous freezing
            !   all the tendencies in those bins are reduced.
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%ms(k,n)%dcondt(j) = gr%ms(k,n)%dcondt(j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            enddo
            enddo
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
            do n = 1, gs%L
            do k = 1, gs%n_bin
              if( icond2(n)==1 ) then
                gs%ms(k,n)%dcondt(j) = gs%ms(k,n)%dcondt(j) * modc_r(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_r(n)
              endif
            end do
            end do
            do n=1,gs%L
              mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
            enddo
          else if( j == 8 ) then
            ! contact ice nucleation
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%ms(k,n)%dcondt(j) = gr%ms(k,n)%dcondt(j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            end do
            end do
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
            do n = 1, gs%L
            do k = 1, gs%n_bin
              if( icond2(n)==1) then
                gs%ms(1,n)%dcondt(j) = gs%ms(1,n)%dcondt(j) * modc_r(n)
                acc_mod_s(1,j,n) = acc_mod_s(1,j,n) * modc_r(n)
              endif
            enddo
            enddo
            do n=1,gs%L
              mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
            enddo
            do n=1,gs%L
              if( icond2(n)==1) then
                ga(2)%ms(1,n)%dcondt(j) = ga(2)%ms(1,n)%dcondt(j) * modc_r(n)
                acc_mod_a(1,j,2,n) = acc_mod_a(1,j,2,n) * modc_r(n)
              endif
              mark_acc_a(2,n)=mark_acc_a(2,n) + icond1(n)*icond2(n)
            enddo
          else if( j == 1 ) then
            ! evaporation
!            do kn=1,gr%n_bin*gr%L
!              n=(kn-1)/gr%n_bin+1
!              k=kn-(n-1)*gr%n_bin
             do n = 1, gr%L
             do k = 1, gr%n_bin
              if( icond2(n)==1) then
                gr%ms(k,n)%dcondt(j) = gr%ms(k,n)%dcondt(j) * modc_r(n)
                acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_r(n)
              endif
            end do
            end do
            do n=1,gr%L
              ick(n)=0
              if(icond2(n)==1.and.gs%MS(i,n)%inevp>0) then
                ick(n)=gr%MS(i,n)%inevp
              endif
            enddo
            do n=1,gr%L
              if(ick(n)>0) then
                ica=ick(n)
                ga(ica)%ms(1,n)%dcondt(j) = ga(ica)%ms(1,n)%dcondt(j) * modc_r(n)
                acc_mod_a(1,j,ica,n) = acc_mod_a(1,j,ica,n) * modc_r(n)
                mark_acc_a(ica,n)=mark_acc_a(ica,n) + icond1(n)
              end if
            enddo
          else
            do n=1,gr%L
              if( icond2(n)==1) then
                gr%MS(i,n)%dcondt(j) = gr%MS(i,n)%dcondt(j) * modc_r(n)
                acc_mod_r(i,j,n) = acc_mod_r(i,j,n) * modc_r(n)
              end if
           end do
          end if
        end do
      enddo liqbin_loop1
    end if
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 3. bugdet of each solid hydrometeor bin
    if( flagp_s > 0 ) then
      rimmlt_if: if(iupdate_gs(4)==1.and.iupdate_gs(10)==1) then

        icebin_loop1: do i=1,gs%n_bin
          do n=1,gs%L
            total_sloss = &
               ! - riming process with rain drops
               max(-gs%MS(i,n)%dcondt(4),0.0_ds) + &
               ! - melting-shedding
               max(-gs%MS(i,n)%dcondt(10),0.0_ds)


            total_sgain = &
               ! - existing con
!!c               0.9999*gs%MS(i,n)%con/gs%dt + &
               gs%MS(i,n)%con/gs%dt+&
               ! - riming process with rain drops
               max(gs%MS(i,n)%dcondt(4),0.0_ds)*real(iupdate_gs(4),PS_KIND)

            icond1(n)=0
!            ierror1(n)=0
            modc_s(n)=1.0_PS

            if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then

              modc_s(n) = min(1.0_ps,&
                  total_sgain/max(total_sloss,total_limit))

              modc_s(n)=max(0.0_PS,modc_s(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_s(n) = mark_acc_s(n) + 1

              icond1(n)=1
!!c             write(*,26) i,n,modc_s,total_sloss
!!c26           format("con:modc_s",2i5,2es15.6)
            endif
            if(modc_s(n)>1.0e+4) then
               LOG_ERROR("cal_con_budget_col",*) "something wrong",n,gs%MS(i,n)%dcondt(1:gs%n_tendpros)
               call PRC_abort
            endif
          enddo

          if(all(icond1(1:gs%L)==0)) cycle

          do j=4,10,6
            do n=1,gs%L
              icond2(n)=0
              if( gs%MS(i,n)%dcondt(j) < 0.0_ps ) then
                icond2(n)=1
              endif
            enddo

            if(j==4) then
              ! riming process with rain drops
              ! ice concentration should not change in total, which means that
              ! if negative tend is decreased, then positive one has to be
              ! decreased.
              !
              do n=1,gs%L
                if(icond2(n)==1) then
                  gs%MS(i,n)%dcondt(j)=gs%MS(i,n)%dcondt(j)*modc_s(n)
                  acc_mod_s(i,j,n)=acc_mod_s(i,j,n)*modc_s(n)
                endif
                total_rims_tend(n) = 0.0_ps
                total_rimr_tend(n) = 0.0_ps
              enddo
              do k=1,gs%n_bin
                do n=1,gs%L
                  total_rims_tend(n)=total_rims_tend(n)+gs%MS(k,n)%dmassdt(imt,j)*acc_mod_s(k,j,n)
                end do
              end do
              do k=1,gr%n_bin
                do n=1,gr%L
                  total_rimr_tend(n)=total_rimr_tend(n)-gr%MS(k,n)%dmassdt(rmt,j)*acc_mod_r(k,j,n)
                end do
              end do
              do n=1,gs%L
                mod_dum(n)=1.0_PS
                if(icond2(n)==1) then
                  mod_dum(n)=max(0.0_PS,min(2.0_PS,total_rims_tend(n)/max(total_limit,total_rimr_tend(n))))
                endif
              enddo

!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                gr%ms(k,n)%dcondt(j)=gr%ms(k,n)%dcondt(j)*mod_dum(n)
                acc_mod_r(k,j,n)=acc_mod_r(k,j,n)*mod_dum(n)
              enddo
              enddo
              do n=1,gr%L
                mark_acc_r(n)=mark_acc_r(n)+icond1(n)*icond2(n)
              enddo
!!c                      write(*,'("ice con>rim",4I5,10ES15.6)') KD(n),ID(n),JD(n),i,&
!!c                             mod_dum,modc_s,total_rims_tend*gs%dt,total_rimr_tend*gs%dt

            elseif(j==10) then
              ! melting-shedding process
              ! total transfer from ice to liq does decrease.
              !
              ! calculate total transfer from ice to liq or liq to ice
              ! all tendencies for melting process is positive for liq and neg for ice.
              do n=1,gs%L
                if(icond2(n)==1) then
                  gs%MS(i,n)%dcondt(j)=gs%MS(i,n)%dcondt(j)*modc_s(n)
                  acc_mod_s(i,j,n)=acc_mod_s(i,j,n)*modc_s(n)
                endif
                total_mlts_tend(n)=0.0_PS
                total_mltr_tend(n)=0.0_PS
              enddo
              do k=1,gs%n_bin
                do n=1,gs%L
                  total_mlts_tend(n)=total_mlts_tend(n)-gs%MS(k,n)%dmassdt(imt,j)*acc_mod_s(k,j,n)
                end do
              end do
              do k=1,gr%n_bin
                do n=1,gs%L
                  total_mltr_tend(n)=total_mltr_tend(n)+gr%MS(k,n)%dmassdt(rmt,j)*acc_mod_r(k,j,n)
                end do
              end do
              do n=1,gs%L
                mod_dum(n)=max(0.0_PS,min(1.0_PS,total_mlts_tend(n)/&
                           max(total_limit,total_mltr_tend(n))))
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%ms(k,n)%dcondt(j)=gr%ms(k,n)%dcondt(j)*mod_dum(n)
                  acc_mod_r(k,j,n)=acc_mod_r(k,j,n)*mod_dum(n)
                endif
              enddo
              enddo
              do n=1,gs%L
                mark_acc_r(n) = mark_acc_r(n) + icond1(n)*icond2(n)
              enddo
            endif
          enddo
        enddo icebin_loop1
      endif rimmlt_if
!!c       total_mlts_tend = 0.0_ps
!!c       total_rims_tend = 0.0_ps
!!c       do k=1,gs%n_bin
!!c          total_mlts_tend=total_mlts_tend-gs%MS(k,n)%dmassdt(imt,10)*acc_mod_s(k,10)
!!c          total_rims_tend=total_rims_tend+gs%MS(k,n)%dmassdt(imt,4)*acc_mod_s(k,4)
!!c       end do
!!c
!!c       total_mltr_tend=0.0
!!c       total_rimr_tend=0.0
!!c       do k=1,gr%n_bin
!!c          total_mltr_tend=total_mltr_tend+gr%MS(k,n)%dmassdt(rmt,10)*acc_mod_r(k,10)
!!c          total_rimr_tend=total_rimr_tend-gr%MS(k,n)%dmassdt(rmt,4)*acc_mod_r(k,4)
!!c       end do
!!c       if(total_mlts_tend>0.0_PS.or.total_rims_tend>0.0_PS) then
!!c          write(*,'("ice con>mltchk",3I5,10ES15.6)') KD(n),ID(n),JD(n),&
!!c                    total_mlts_tend*gr%dt,total_mltr_tend*gr%dt,&
!!c                    total_rims_tend*gr%dt,total_rimr_tend*gr%dt
!!c       endif

      icebin_loop2: do i=1,gs%n_bin
        do n=1,gs%L
          total_sloss = &
               ! - collision advection
               max(-gs%MS(i,n)%dcondt(2),0.0_ds)*real(iupdate_gs(2),PS_KIND) + &
               ! - riming process with cloud droplets
               max(-gs%MS(i,n)%dcondt(3),0.0_ds)*real(iupdate_gs(3),PS_KIND) + &
               ! - riming process with rain drops
               max(-gs%MS(i,n)%dcondt(4),0.0_ds)*real(iupdate_gs(4),PS_KIND) + &
               ! - collision breakup
               max(-gs%MS(i,n)%dcondt(5),0.0_ds)*real(iupdate_gs(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(-gs%MS(i,n)%dcondt(6),0.0_ds)*real(iupdate_gs(6),PS_KIND) + &
               ! - secondary nucleation (splintering)
               max(-gs%MS(i,n)%dcondt(9), 0.0_ds)*real(iupdate_gs(9),PS_KIND) + &
               ! - melting-shedding
               max(-gs%MS(i,n)%dcondt(10),0.0_ds)*real(iupdate_gs(10),PS_KIND)


          total_sgain = &
               ! - existing con
!!c               0.9999*gs%MS(i,n)%con/gs%dt + &
               gs%MS(i,n)%con/gs%dt + &
               ! - collision advection
               max(gs%MS(i,n)%dcondt(2),0.0_ds)*real(iupdate_gs(2),PS_KIND) + &
               ! - riming process with cloud droplets
               max(gs%MS(i,n)%dcondt(3),0.0_ds)*real(iupdate_gs(3),PS_KIND) + &
               ! - riming process with rain drops
               max(gs%MS(i,n)%dcondt(4),0.0_ds)*real(iupdate_gs(4),PS_KIND) + &
               ! - collision breakup
               max(gs%MS(i,n)%dcondt(5),0.0_ds)*real(iupdate_gs(5),PS_KIND) + &
               ! - hydrodynamic breakup
               max(gs%MS(i,n)%dcondt(6),0.0_ds)*real(iupdate_gs(6),PS_KIND) + &
               ! - ice nucleation (contact)
               max(gs%MS(i,n)%dcondt(8),0.0_ds)*real(iupdate_gs(8),PS_KIND) + &
               ! - secondary nucleation (splintering)
               max( gs%MS(i,n)%dcondt(9), 0.0_ds)*real(iupdate_gs(9),PS_KIND) + &
               ! - ice nucleation (immersion)
               max( gs%MS(i,n)%dcondt(11), 0.0_ds)*real(iupdate_gs(11),PS_KIND) + &
               ! - ice nucleation (homogeneous freezing)
               max( gs%MS(i,n)%dcondt(12), 0.0_ds)*real(iupdate_gs(12),PS_KIND)

!          ierror1(n)=0
          icond1(n)=0
          modc_s(n)=1.0_PS

          if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then

            dum1=max(-gs%MS(i,n)%dcondt(4),0.0_DS)*real(iupdate_gs(4),PS_KIND)&
                   +max(-gs%MS(i,n)%dcondt(10),0.0_DS)*real(iupdate_gs(10),PS_KIND)

            modc_s(n) = min(1.0_ps,&
                (total_sgain-dum1)/max(total_sloss-dum1,total_limit))

            if(total_sloss-dum1<total_limit) then
              modc_s(n)=1.0_PS
            else
              modc_s(n)=max(0.0_PS,modc_s(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_s(n) = mark_acc_s(n) + 1
              icond1(n)=1
            endif
          else if( total_sgain == 0.0_ps .and. total_sloss == 0.0_ps ) then
            ! this mark should help having bogus tendency in
            ! calculating new tendency (cal_model_tendency)
            if( gs%MS(i,n)%mark == 4 ) then
              gs%MS(i,n)%mark = 6
            else
              gs%MS(i,n)%mark = 5
            end if
          end if
          if(modc_s(n)>1.0e+4) then
             LOG_ERROR("cal_con_budget_col",*) "something wrong",n,gs%MS(i,n)%dcondt(1:gs%n_tendpros)
              call PRC_abort
           endif
        enddo

        if(all(icond1(1:gs%L)==0)) cycle

        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
          do n=1,gs%L
            icond2(n)=0
            if( gs%MS(i,n)%dcondt(j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo
          if(j == 2 .or. j == 5 .or. j == 6 .or. j==9) then
            ! - if those processes contribute to the total loss,
            !   the tendencies in all the bins are reduced for
            !   consistency.
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
             do n = 1, gs%L
             do k = 1, gs%n_bin
              if(icond2(n)==1) then
                gs%ms(k,n)%dcondt(j) = gs%ms(k,n)%dcondt(j) * modc_s(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_s(n)
              endif
            end do
            end do
          else if( j == 1 ) then
            ! evaporation
             LOG_ERROR("cal_con_budget_col",*) "wrong in con_col",j
             call PRC_abort
          else
            do n=1,gs%L
              if(icond2(n)==1) then
                gs%MS(i,n)%dcondt(j) = gs%MS(i,n)%dcondt(j) * modc_s(n)
                acc_mod_s(i,j,n) = acc_mod_s(i,j,n) * modc_s(n)
              end if
            end do
          endif
        end do
      end do icebin_loop2
    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 4. bugdet of each aerosol particles
    if( flagp_a /= 0 ) then
      do ica=1,ncat_a
        aerbin_loop1: do i=1,ga(ica)%n_bin
          do n=1,ga(ica)%L
            total_aloss = &
               ! - ccn activation -
               max(-ga(ica)%MS(i,n)%dcondt(3),0.0_ds)*real(iupdate_ga(3,ica),PS_KIND) + &
               ! - ice nucleation (deposition/sorption)
               max(-ga(ica)%MS(i,n)%dcondt(7),0.0_ds)*real(iupdate_ga(7,ica),PS_KIND) + &
               ! - ice nucleation (contact)
               max(-ga(ica)%MS(i,n)%dcondt(8),0.0_ds)*real(iupdate_ga(8,ica),PS_KIND)

            exist_con=max(ga(ica)%MS(i,n)%con-n_lmt_ap,0.0_ps)

            total_again = &
               ! - existing con
!!c                  0.9999*exist_con/ga(ica)%dt + &
               exist_con/ga(ica)%dt + &
               ! - vapor deposition (evaporation)
               max(ga(ica)%MS(i,n)%dcondt(1),0.0_ds)*real(iupdate_ga(1,ica),PS_KIND)

!            ierror1(n)=0
            icond1(n)=0
            modc_a(n)=1.0_PS

            if( total_again < 0.99999_ps*total_aloss.and.total_aloss>total_limit ) then

              dum1=max(-ga(ica)%MS(i,n)%dcondt(3),0.0_DS)*real(iupdate_ga(3,ica),PS_KIND) &
                  +max(-ga(ica)%MS(i,n)%dcondt(7),0.0_DS)*real(iupdate_ga(7,ica),PS_KIND)

              modc_a(n) = min(1.0_ps,&
                    (total_again-dum1)/max(total_aloss-dum1,total_limit))

              if(total_aloss-dum1<total_limit) then
                modc_a(n)=1.0_PS
              else
                modc_a(n)=max(0.0_PS,modc_a(n))

                mark_mod(n) = mark_mod(n) + 1
                mark_acc_a(ica,n) = mark_acc_a(ica,n) + 1

                icond1(n)=1
              endif
            else if( total_again == 0.0_ps .and. total_aloss == 0.0_ps ) then
              ! this mark should help having bogus tendency in
              ! calculating new tendency (cal_model_tendency)
              if( ga(ica)%MS(i,n)%mark == 4 ) then
                ga(ica)%MS(i,n)%mark = 6
              else
                ga(ica)%MS(i,n)%mark = 5
              end if
            end if
            if(modc_a(n)>1.0e+4_ps) then
               LOG_ERROR("cal_con_budget_col",*) "something wrong",n,ga(ica)%MS(i,n)%dcondt(1:ga(ica)%n_tendpros)
                call PRC_abort
            endif
            if( (total_aloss==&
                    max(-ga(ica)%MS(i,n)%dcondt(3),0.0_DS)*real(iupdate_ga(3,ica),PS_KIND) &
                  + max(-ga(ica)%MS(i,n)%dcondt(7),0.0_DS)*real(iupdate_ga(7,ica),PS_KIND) ) &
                  .and. icond1(n) > 0 ) then
               LOG_ERROR("cal_con_budget_col",*) "something is not right at cal_con_budet,i,n",i,n
               LOG_ERROR_CONT(*) "ica,total_aloss,again",ica,total_aloss,total_again
               LOG_ERROR_CONT(*) icond1(n)
                do j = 1, ga(ica)%n_tendpros
                   LOG_ERROR_CONT(*) j,ga(ica)%MS(i,n)%dcondt(j)
                end do
                call PRC_abort
            endif
          enddo

          if(all(icond1(1:ga(ica)%L)==0)) cycle

          do j=1,ga(ica)%n_tendpros
            if(iupdate_ga(j,ica)==0) cycle
            do n=1,ga(ica)%L
              icond2(n)=0
              if( ga(ica)%MS(i,n)%dcondt(j) < 0.0_ps ) then
                icond2(n)=1
              endif
            enddo

            if(j==3) then
              ! ccn activation
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dcondt(j) = ga(ica)%ms(k,n)%dcondt(j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%ms(k,n)%dcondt(j) = gr%ms(k,n)%dcondt(j) * modc_a(n)
                  acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gr%L
                mark_acc_r(n) = mark_acc_r(n) + icond1(n)*icond2(n)
              enddo

            elseif(j==7) then
              ! ice nucleation (deposition/sorption)
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dcondt(j) = ga(ica)%ms(k,n)%dcondt(j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              enddo
              enddo
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%ms(k,n)%dcondt(j) = gs%ms(k,n)%dcondt(j) * modc_a(n)
                  acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gs%L
                mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
              enddo
            elseif(j==8) then
              ! ice nucleation (contact)
!              do kn=1,ga(ica)%n_bin*ga(ica)%L
!                n=(kn-1)/ga(ica)%n_bin+1
!                k=kn-(n-1)*ga(ica)%n_bin
               do n = 1, ga(ica)%L
               do k = 1, ga(ica)%n_bin
                if(icond2(n)==1) then
                  ga(ica)%ms(k,n)%dcondt(j) = ga(ica)%ms(k,n)%dcondt(j) * modc_a(n)
                  acc_mod_a(k,j,ica,n) = acc_mod_a(k,j,ica,n) * modc_a(n)
                endif
              end do
              end do
!              do kn=1,gs%n_bin*gs%L
!                n=(kn-1)/gs%n_bin+1
!                k=kn-(n-1)*gs%n_bin
              do n = 1, gs%L
              do k = 1, gs%n_bin
                if(icond2(n)==1) then
                  gs%ms(k,n)%dcondt(j) = gs%ms(k,n)%dcondt(j) * modc_a(n)
                  acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_a(n)
                endif
              end do
              end do
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                if(icond2(n)==1) then
                  gr%MS(k,n)%dmassdt(rmt,j)=gr%MS(k,n)%dmassdt(rmt,j) * modc_a(n)
                  acc_mod_r(k,j,n) = acc_mod_r(k,j,n) * modc_a(n)
                endif
              enddo
              enddo
              do n=1,gr%L
                mark_acc_r(n) = mark_acc_r(n) + icond1(n)*icond2(n)
                mark_acc_s(n) = mark_acc_s(n) + icond1(n)*icond2(n)
              enddo
            endif
          enddo
        enddo aerbin_loop1
      enddo
    endif
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine cal_con_budget_col

  subroutine cal_vol_budget_col( gr, gs, nvol,&
       acc_mod_r, acc_mod_s,acc_mod_a, &
       mark_mod, mark_acc_r , mark_acc_s, mark_acc_a,flagp_s,&
       iupdate_gs&
       )
    use scale_prc, only: &
       PRC_abort
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate total gain and loss for each group, and modify
    ! only concentration tendency
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (group), intent(inout)         :: gr, gs
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_r
    real (PS), dimension(mxnbin,mxntend,*)    :: acc_mod_s
    real (PS), dimension(mxnbina,mxntend,ncamx,*)    :: acc_mod_a
    integer, dimension(*), intent(inout)     :: mark_mod
    integer, dimension(*),intent(inout)     :: mark_acc_r, mark_acc_s
    integer, dimension(ncamx,*)    :: mark_acc_a
    integer, intent(in)                   :: flagp_s!, flagp_a, flagp_r
    ! argument position of volume tendency
    integer, intent(in)                 :: nvol
    !integer, intent(in)                 :: ncat_a
    !integer :: ID(*),JD(*),KD(*)
    integer,dimension(*),intent(in) :: iupdate_gs!,iupdate_gr
    !
    ! local space
    !
    real(ps),dimension(LMAX)    :: total_rimr_tend, total_rims_tend,total_mlts_tend,total_mltr_tend
    real(ps)    :: total_sloss, total_sgain
    ! modification constant for solid_hydro
    real(ps),dimension(LMAX)    :: modc_s
    real(ps)    :: dum1
    real(ps),dimension(LMAX)    :: mod_dum
    integer,dimension(LMAX) :: ick
    integer,dimension(LMAX) :: icond1,icond2 !,ierror1
    ! existing volumes
    real(ps)    :: exist_vol
    real(ps),parameter :: total_limit=1.0e-30_ps
    integer      :: i,j,k,n,ica!,in

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize
    do n=1,gr%L
      mark_mod(n) = 0
    enddo
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 3. bugdet of each solid hydrometeor bin
    if( flagp_s > 0 ) then
      rimmlt_if: if(iupdate_gs(4)==1.and.iupdate_gs(10)==1) then
        icebin_loop1: do i=1,gs%n_bin
          do n=1,gs%L

            total_sloss = &
               ! - riming process with rain drops
               max(-gs%MS(i,n)%dvoldt(nvol,4),0.0_ds) + &
               ! - melting-shedding
               max(-gs%MS(i,n)%dvoldt(nvol,10),0.0_ds)

            select case(nvol)
            case(2)
              exist_vol = gs%MS(i,n)%con*gs%MS(i,n)%a_len**3.0
            case(3)
              exist_vol = gs%MS(i,n)%con*gs%MS(i,n)%c_len**3.0
            case default
              exist_vol = gs%MS(i,n)%con*gs%IS(i,n)%v_cs
            end select


            total_sgain = &
               ! - existing volume
!!c               0.9999*exist_vol/gs%dt + &
               exist_vol/gs%dt + &
               ! - riming process with rain drops
               max(gs%MS(i,n)%dvoldt(nvol,4),0.0_ds) + &
               ! - melting-shedding
               max(gs%MS(i,n)%dvoldt(nvol,10),0.0_ds)

!            ierror1(n)=0
            icond1(n)=0
            modc_s(n)=1.0_PS

            if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then


              modc_s(n) = min(1.0_ps,&
                total_sgain/max(total_sloss,total_limit))

              modc_s(n)=max(0.0_PS,modc_s(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_s(n) = mark_acc_s(n) + 1
              icond1(n)=1

            endif
            if(modc_s(n)>1.0e+4) then
               LOG_ERROR("cal_vol_budget_col",*) "something wrong",gs%MS(i,n)%dvoldt(nvol,1:gs%n_tendpros)
               call PRC_abort
            endif
          enddo

          if(all(icond1(1:gs%L)==0)) cycle

          do j=4,10,6
            do n=1,gs%L
              icond2(n)=0
              if( gs%MS(i,n)%dvoldt(nvol,j) < 0.0_ps ) then
                icond2(n)=1
              endif
            enddo

            if(j==4) then
              ! riming process with rain drops
              !
              do n=1,gs%L
                if(icond2(n)==1) then
                  gs%MS(i,n)%dvoldt(nvol,j)=gs%MS(i,n)%dvoldt(nvol,j)*modc_s(n)
                  acc_mod_s(i,j,n)=acc_mod_s(i,j,n)*modc_s(n)
                endif
              enddo
              do n=1,gs%L
                total_rims_tend(n) = 0.0_ps
                total_rimr_tend(n) = 0.0_ps
              enddo
              do k=1,gs%n_bin
                do n=1,gs%L
                  total_rims_tend(n)=total_rims_tend(n)+gs%MS(k,n)%dmassdt(imt,j)*acc_mod_s(k,j,n)
                enddo
              enddo
              do k=1,gr%n_bin
                do n=1,gr%L
                  total_rimr_tend(n)=total_rimr_tend(n)-gr%MS(k,n)%dmassdt(rmt,j)*acc_mod_r(k,j,n)
                enddo
              enddo
              do n=1,gs%L
                mod_dum(n)=1.0_PS
                if(icond2(n)==1) then
                  mod_dum(n)=max(0.0_PS,min(1.0_PS,total_rims_tend(n)/max(total_limit,total_rimr_tend(n))))
                endif
              enddo

!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                acc_mod_r(k,j,n)=acc_mod_r(k,j,n)*mod_dum(n)
              end do
              end do
              do n=1,gr%L
                mark_acc_r(n)=mark_acc_r(n)+icond1(n)*icond2(n)
              enddo

            elseif(j==10) then
              ! melting-shedding process
              ! total transfer from ice to liq does decrease.
              !
              ! calculate total transfer from ice to liq or liq to ice
              ! all tendencies for melting process is positive for liq and neg for ice.
              do n=1,gs%L
                if(icond2(n)==1) then
                  gs%MS(i,n)%dvoldt(nvol,j)=gs%MS(i,n)%dvoldt(nvol,j)*modc_s(n)
                  acc_mod_s(i,j,n)=acc_mod_s(i,j,n)*modc_s(n)
                endif
                total_mlts_tend(n) = 0.0_ps
                total_mltr_tend(n) = 0.0_ps
              enddo
              do k=1,gs%n_bin
                do n=1,gs%L
                  total_mlts_tend(n)=total_mlts_tend(n)-gs%MS(k,n)%dmassdt(imt,j)*acc_mod_s(k,j,n)
                enddo
              enddo
              do k=1,gr%n_bin
                do n=1,gr%L
                  total_mltr_tend(n)=total_mltr_tend(n)+gr%MS(k,n)%dmassdt(rmt,j)*acc_mod_r(k,j,n)
                enddo
              enddo

              do n=1,gs%L
                mod_dum(n)=1.0_PS
                if(icond2(n)==1) then
                  mod_dum(n)=max(0.0_PS,min(1.0_PS,total_mlts_tend(n)/max(total_limit,total_mltr_tend(n))))
                endif
              enddo
!              do kn=1,gr%n_bin*gr%L
!                n=(kn-1)/gr%n_bin+1
!                k=kn-(n-1)*gr%n_bin
              do n = 1, gr%L
              do k = 1, gr%n_bin
                acc_mod_r(k,j,n)=acc_mod_r(k,j,n)*mod_dum(n)
              end do
              end do
              do n=1,gr%L
                mark_acc_r(n)=mark_acc_r(n)+icond1(n)*icond2(n)
              enddo
            endif
          enddo
        enddo icebin_loop1
      endif rimmlt_if

      icebin_loop2: do i=1,gs%n_bin
        do n=1,gs%L
          total_sloss = &
            ! - collision advection
            max(-gs%MS(i,n)%dvoldt(nvol,2),0.0_ds)*real(iupdate_gs(2),PS_KIND) + &
            ! - riming process with cloud droplets
            max(-gs%MS(i,n)%dvoldt(nvol,3),0.0_ds)*real(iupdate_gs(3),PS_KIND) + &
            ! - riming process with rain drops
            max(-gs%MS(i,n)%dvoldt(nvol,4),0.0_ds)*real(iupdate_gs(4),PS_KIND) + &
            ! - collision breakup
            max(-gs%MS(i,n)%dvoldt(nvol,5),0.0_ds)*real(iupdate_gs(5),PS_KIND) + &
            ! - hydrodynamic breakup
            max(-gs%MS(i,n)%dvoldt(nvol,6),0.0_ds)*real(iupdate_gs(6),PS_KIND) + &
            ! - secondary nucleation (splintering)
            max(-gs%MS(i,n)%dvoldt(nvol,9), 0.0_ds)*real(iupdate_gs(9),PS_KIND) + &
            ! - melting-shedding
            max(-gs%MS(i,n)%dvoldt(nvol,10),0.0_ds)*real(iupdate_gs(10),PS_KIND)

          select case(nvol)
          case(2)
            exist_vol = gs%MS(i,n)%con*gs%MS(i,n)%a_len**3.0
          case(3)
            exist_vol = gs%MS(i,n)%con*gs%MS(i,n)%c_len**3.0
          case default
            exist_vol = gs%MS(i,n)%con*gs%IS(i,n)%v_cs
          end select


          total_sgain = &
          ! - existing volume
            exist_vol/gs%dt + &
            ! - collision advection
            max(gs%MS(i,n)%dvoldt(nvol,2),0.0_ds)*real(iupdate_gs(2),PS_KIND) + &
            ! - riming process with cloud droplets
            max(gs%MS(i,n)%dvoldt(nvol,3),0.0_ds)*real(iupdate_gs(3),PS_KIND) + &
            ! - riming process with rain drops
            max(gs%MS(i,n)%dvoldt(nvol,4),0.0_ds)*real(iupdate_gs(4),PS_KIND) + &
            ! - collision breakup
            max(gs%MS(i,n)%dvoldt(nvol,5),0.0_ds)*real(iupdate_gs(5),PS_KIND) + &
            ! - hydrodynamic breakup
            max(gs%MS(i,n)%dvoldt(nvol,6),0.0_ds)*real(iupdate_gs(6),PS_KIND) + &
            ! - ice nucleation (contact)
            max(gs%MS(i,n)%dvoldt(nvol,8),0.0_ds)*real(iupdate_gs(8),PS_KIND) + &
            ! - secondary nucleation (splintering)
            max(gs%MS(i,n)%dvoldt(nvol,9), 0.0_ds)*real(iupdate_gs(9),PS_KIND) + &
            ! - melting-shedding
            max(gs%MS(i,n)%dvoldt(nvol,10),0.0_ds)*real(iupdate_gs(10),PS_KIND) + &
            ! - ice nucleation (immersion)
            max(gs%MS(i,n)%dvoldt(nvol,11), 0.0_ds)*real(iupdate_gs(11),PS_KIND) + &
            ! - ice nucleation (homogeneous freezing)
            max(gs%MS(i,n)%dvoldt(nvol,12), 0.0_ds)*real(iupdate_gs(12),PS_KIND)

!          ierror1(n)=0
          icond1(n)=0
          modc_s(n)=1.0_PS

          if( total_sgain < 0.99999_ps*total_sloss.and.total_sloss>total_limit ) then

            dum1=max(-gs%MS(i,n)%dvoldt(nvol,4),0.0_DS)*real(iupdate_gs(4),PS_KIND) &
                +max(-gs%MS(i,n)%dvoldt(nvol,10),0.0_DS)*real(iupdate_gs(10),PS_KIND)

            modc_s(n) = min(1.0_ps,&
                (total_sgain-dum1)/max(total_sloss-dum1,total_limit))

            if(total_sloss-dum1<total_limit) then
              modc_s(n)=1.0_PS
            else
              modc_s(n)=max(0.0_PS,modc_s(n))

              mark_mod(n) = mark_mod(n) + 1
              mark_acc_s(n) = mark_acc_s(n) + 1
              icond1(n)=1
            endif
          end if
          if(modc_s(n)>1.0e+4) then
             LOG_ERROR("cal_vol_budget_col",*) "something wrong",i,n,nvol,gs%MS(i,n)%dvoldt(nvol,1:gs%n_tendpros)
             call PRC_abort
           endif
        enddo
!!c             write(*,27) KD(n),ID(n),JD(n),i,modc_s,total_sloss,total_sgain,&
!!c                        gs%MS(i,n)%dvoldt(nvol,1),gs%MS(i,n)%dvoldt(nvol,10)
!!c27        format("vol:modc_s 2",4i5,20es15.6)

        if(all(icond1(1:gs%L)==0)) cycle

        do j=1,gs%n_tendpros
          if(iupdate_gs(j)==0) cycle
          do n=1,gs%L
            icond2(n)=0
            if( gs%MS(i,n)%dvoldt(nvol,j) < 0.0_ps ) then
              icond2(n)=1
            endif
          enddo

          if(j == 2 .or. j == 5 .or. j == 6 .or. j==9) then
            ! - if those processes contribute to the total loss,
            !   the tendencies in all the bins are reduced for
            !   consistency.
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
             do n = 1, gs%L
             do k = 1, gs%n_bin
              if(icond2(n)==1) then
                gs%ms(k,n)%dvoldt(nvol,j) = gs%ms(k,n)%dvoldt(nvol,j) * modc_s(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_s(n)
              endif
            end do
            end do
          else if( j == 1 ) then
            ! evaporation
!            do kn=1,gs%n_bin*gs%L
!              n=(kn-1)/gs%n_bin+1
!              k=kn-(n-1)*gs%n_bin
             do n = 1, gs%L
             do k = 1, gs%n_bin
              if(icond2(n)==1) then
                gs%ms(k,n)%dvoldt(nvol,j) = gs%ms(k,n)%dvoldt(nvol,j) * modc_s(n)
                acc_mod_s(k,j,n) = acc_mod_s(k,j,n) * modc_s(n)
              endif
            end do
            end do
            do n=1,gs%L
              ick(n)=0
              if(icond2(n)==1.and.gs%MS(i,n)%inevp>0) then
                ick(n)=gs%MS(i,n)%inevp
              endif
            enddo
            do n=1,gs%L
              if(ick(n)>0) then
                ica=ick(n)
                acc_mod_a(1,j,ica,n) = acc_mod_a(1,j,ica,n) * modc_s(n)
                mark_acc_a(ica,n)=mark_acc_a(ica,n) + icond1(n)
              end if
            enddo
          else
            do n=1,gs%L
              if( icond2(n)==1) then
                gs%MS(i,n)%dvoldt(nvol,j) = gs%MS(i,n)%dvoldt(nvol,j) * modc_s(n)
                acc_mod_s(i,j,n) = acc_mod_s(i,j,n) * modc_s(n)
              endif
            enddo
          endif
        enddo
      enddo icebin_loop2

    endif
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end subroutine cal_vol_budget_col


  subroutine mod_other_tendency_vap(level, g, acc_mod, switch)
    use scale_prc, only: &
       PRC_abort
    integer,intent(in) :: level
    type (group), intent(inout)         :: g
    real (ps), dimension(mxnbin,mxntend,*)    :: acc_mod
    ! switch to calculate the tendency
    ! 1: tendencies of concentration and volume are calculated.
    ! 2: tendencies of mass and volume are calculated.
    integer   :: switch
    integer   :: i,j,k,n

    select case(switch)
    case (1)
      if(g%token==1) then
        do j=1,3,2
          ! vapor deposition (1) and CCN activation (3) only
          do k = 1, g%n_masscom
            if(level>=4.and.(j==1.or.j==3).and.k==rmat-1) cycle
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
            do n = 1, g%L
            do i = 1, g%N_BIN
              g%MS(i,n)%dmassdt(1+k,j) = g%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,n)
            end do
            end do
          enddo

!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
          do n = 1, g%L
          do i = 1, g%N_BIN
            g%MS(i,n)%dcondt(j) = g%MS(i,n)%dcondt(j)*acc_mod(i,j,n)
          enddo
          enddo

!          do k = 1, g%n_nonmass
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
          do n = 1, g%L
          do i = 1, g%N_BIN
          do k = 1, g%n_nonmass
              g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
          enddo
          enddo
          enddo
        enddo

      elseif(g%token==2) then
        do j=1,g%n_tendpros
          if(j/=1.or.j/=3.or.j/=7) cycle

          ! vapor deposition (1), IN (dep) (7) and IN (dhf) (3)
!          do k = 1, g%n_masscom
!            if(level>=4.and.(j==1.or.j==3.or.j==7).and.k==imat-1) cycle
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
          do n = 1, g%L
          do i = 1, g%N_BIN
             do k = 1, g%n_masscom
                if(level>=4.and.(j==1.or.j==3.or.j==7).and.k==imat-1) cycle
                g%MS(i,n)%dmassdt(1+k,j) = g%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,n)
             end do
          enddo
          enddo

!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
          do n = 1, g%L
          do i = 1, g%N_BIN
             g%MS(i,n)%dcondt(j) = g%MS(i,n)%dcondt(j)*acc_mod(i,j,n)
          enddo
          enddo

!          do k = 1, g%n_nonmass
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
          do n = 1, g%L
          do i = 1, g%N_BIN
          do k = 1, g%n_nonmass
             g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
          enddo
          enddo
          enddo
        enddo
     end if

    case default
       LOG_ERROR("mod_other_tendency_vap",*) "Not ready yet",switch
       call PRC_abort
    end select

  end subroutine mod_other_tendency_vap

  subroutine mod_other_tendency(g, acc_mod, switch, iupdate)
    use scale_prc, only: &
       PRC_abort
    !integer,intent(in) :: level
    type (group), intent(inout)         :: g
    real (ps), dimension(mxnbin,mxntend,*)    :: acc_mod
    ! switch to calculate the tendency
    ! 1: tendencies of concentration and volume are calculated.
    ! 2: tendencies of mass and volume are calculated.
    integer,intent(in)   :: switch
    integer,dimension(*),intent(in) :: iupdate
    integer   :: i, j, k,n

    select case(switch)
    case (1)
      do j = 1, g%n_tendpros
        if(iupdate(j)==0) cycle

!        do k = 1, g%n_masscom
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, g%n_masscom
           g%MS(i,n)%dmassdt(1+k,j) = g%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,n)
        enddo
        enddo
        enddo
!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
          g%MS(i,n)%dcondt(j) = g%MS(i,n)%dcondt(j)*acc_mod(i,j,n)
        enddo
        enddo
!        do k = 1, g%n_nonmass
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, g%n_nonmass
           g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
        end do
        end do
        end do
      end do
    case (2)
      do j = 1, g%n_tendpros
        if(iupdate(j)==0) cycle
!        do k = 1, 1+g%n_masscom
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, 1+g%n_masscom
           g%MS(i,n)%dmassdt(k,j) = g%MS(i,n)%dmassdt(k,j)*acc_mod(i,j,n)
        end do
        end do
        end do

!        do k = 1, g%n_nonmass
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, g%n_nonmass
           g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
        end do
        end do
        end do
      end do
    case (3)
      do j = 1, g%n_tendpros
        if(iupdate(j)==0) cycle
!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
          g%MS(i,n)%dcondt(j) = g%MS(i,n)%dcondt(j)*acc_mod(i,j,n)
        enddo
        enddo
!        do k = 1, 1+g%n_masscom
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, 1+g%n_masscom
           g%MS(i,n)%dmassdt(k,j) = g%MS(i,n)%dmassdt(k,j)*acc_mod(i,j,n)
        end do
        end do
        end do
        if(g%token==2) then
!          do k = 1, g%n_nonmass
!            if( k /= ivcs ) then
!              do in=1,g%n_bin*g%L
!                n=(in-1)/g%N_BIN+1
!                i=in-(n-1)*g%N_BIN
           do n = 1, g%L
           do i = 1, g%N_BIN
           do k = 1, g%n_nonmass
              if( k /= ivcs ) then
                 g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
              endif
           end do
           end do
           end do
        else
!          do k = 1, g%n_nonmass
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
           do n = 1, g%L
           do i = 1, g%N_BIN
           do k = 1, g%n_nonmass
              g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
           end do
           end do
           end do
        end if
      end do

    case (4)
      do j = 1, g%n_tendpros
        if(iupdate(j)==0) cycle
!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
          g%MS(i,n)%dcondt(j) = g%MS(i,n)%dcondt(j)*acc_mod(i,j,n)
        enddo
        enddo
!        do k = 1, 1+g%n_masscom
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, 1+g%n_masscom
           g%MS(i,n)%dmassdt(k,j) = g%MS(i,n)%dmassdt(k,j)*acc_mod(i,j,n)
        end do
        end do
        end do
        if(g%token==2) then
!          do k = 1, g%n_nonmass
!            if( k /= iacr ) then
!              do in=1,g%n_bin*g%L
!                n=(in-1)/g%N_BIN+1
!                i=in-(n-1)*g%N_BIN
           do n = 1, g%L
           do i = 1, g%N_BIN
           do k = 1, g%n_nonmass
              if( k /= iacr ) then
                 g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
              end if
           end do
           end do
           end do
        else
!          do k = 1, g%n_nonmass
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
           do n = 1, g%L
           do i = 1, g%N_BIN
           do k = 1, g%n_nonmass
              g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
           end do
           end do
           end do
        endif
      end do
    case (5)
      do j = 1, g%n_tendpros
        if(iupdate(j)==0) cycle
!        do in=1,g%n_bin*g%L
!          n=(in-1)/g%N_BIN+1
!          i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
          g%MS(i,n)%dcondt(j) = g%MS(i,n)%dcondt(j)*acc_mod(i,j,n)
        enddo
        enddo
!        do k = 1, 1+g%n_masscom
!          do in=1,g%n_bin*g%L
!            n=(in-1)/g%N_BIN+1
!            i=in-(n-1)*g%N_BIN
        do n = 1, g%L
        do i = 1, g%N_BIN
        do k = 1, 1+g%n_masscom
           g%MS(i,n)%dmassdt(k,j) = g%MS(i,n)%dmassdt(k,j)*acc_mod(i,j,n)
        end do
        end do
        end do
        if(g%token==2) then
!          do k = 1, g%n_nonmass
!            if( k /= iccr ) then
!              do in=1,g%n_bin*g%L
!                n=(in-1)/g%N_BIN+1
!                i=in-(n-1)*g%N_BIN
           do n = 1, g%L
           do i = 1, g%N_BIN
           do k = 1, g%n_nonmass
              if( k /= iccr ) then
                 g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
              end if
           end do
           end do
           end do
        else
!          do k = 1, g%n_nonmass
!            do in=1,g%n_bin*g%L
!              n=(in-1)/g%N_BIN+1
!              i=in-(n-1)*g%N_BIN
           do n = 1, g%L
           do i = 1, g%N_BIN
           do k = 1, g%n_nonmass
              g%MS(i,n)%dvoldt(k,j) = g%MS(i,n)%dvoldt(k,j)*acc_mod(i,j,n)
           end do
           end do
           end do
        endif
      end do
    case default
       LOG_ERROR("mod_other_tendency",*) "Not ready yet",switch
       call PRC_abort
    end select

  end subroutine mod_other_tendency

  subroutine mod_other_tendency_ap( g, acc_mod, ica, switch, iupdate_ga)
    use scale_prc, only: &
       PRC_abort
    type (group), dimension(*)         :: g
!tmp    real (ps), pointer, dimension(:,:,:)    :: acc_mod
    integer,intent(in) :: ica!,ncat_a
    real (ps), dimension(mxnbina,mxntend,ncamx,*)    :: acc_mod
    ! switch to calculate the tendency
    ! 1: tendencies of concentration and volume are calculated.
    ! 2: tendencies of mass and volume are calculated.
    integer,intent(in)   :: switch
    integer,dimension(mxntend,*),intent(in) :: iupdate_ga
    integer   :: i, j, k,n

    select case(switch)
    case (1)
      do j=1,g(ica)%n_tendpros
        if(iupdate_ga(j,ica)==0) cycle
!        do k = 1, g(ica)%n_masscom
!          do in=1,g(ica)%n_bin*g(ica)%L
!            n=(in-1)/g(ica)%N_BIN+1
!            i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
        do k = 1, g(ica)%n_masscom
           g(ica)%MS(i,n)%dmassdt(1+k,j) = g(ica)%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,ica,n)
        end do
        end do
        end do
!        do in=1,g(ica)%n_bin*g(ica)%L
!          n=(in-1)/g(ica)%N_BIN+1
!          i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
          g(ica)%MS(i,n)%dcondt(j) = g(ica)%MS(i,n)%dcondt(j)*acc_mod(i,j,ica,n)
        end do
        end do
      end do
    case (2)
      do j=1,g(ica)%n_tendpros
        if(iupdate_ga(j,ica)==0) cycle
!        do k = 1, 1+g(ica)%n_masscom
!          do in=1,g(ica)%n_bin*g(ica)%L
!            n=(in-1)/g(ica)%N_BIN+1
!            i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
        do k = 1, 1+g(ica)%n_masscom
           g(ica)%MS(i,n)%dmassdt(k,j) = g(ica)%MS(i,n)%dmassdt(k,j)*acc_mod(i,j,ica,n)
        end do
        end do
        end do
      end do
    case default
       LOG_ERROR("mod_other_tendency_ap",*) "Not ready yet",switch
       call PRC_abort
    end select

  end subroutine mod_other_tendency_ap

  subroutine mod_other_tendency_ap_vap( g, acc_mod, ica, switch)
    use scale_prc, only: &
       PRC_abort
    type (group), dimension(*)         :: g
    integer,intent(in) :: ica!,ncat_a
    real (ps), dimension(mxnbina,mxntend,ncamx,*)    :: acc_mod
    ! switch to calculate the tendency
    ! 1: tendencies of concentration and volume are calculated.
    ! 2: tendencies of mass and volume are calculated.
    integer   :: switch
    integer   :: i, j, k, n

    select case(switch)
    case (1)
      j=1
!      do k = 1, g(ica)%n_masscom
!        do in=1,g(ica)%n_bin*g(ica)%L
!          n=(in-1)/g(ica)%N_BIN+1
!          i=in-(n-1)*g(ica)%N_BIN
      do n = 1, g(ica)%L
      do i = 1, g(ica)%N_BIN
      do k = 1, g(ica)%n_masscom
         g(ica)%MS(i,n)%dmassdt(1+k,j) = g(ica)%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,ica,n)
      end do
      end do
      end do
!      do in=1,g(ica)%n_bin*g(ica)%L
!        n=(in-1)/g(ica)%N_BIN+1
!        i=in-(n-1)*g(ica)%N_BIN
      do n = 1, g(ica)%L
      do i = 1, g(ica)%N_BIN
        g(ica)%MS(i,n)%dcondt(j) = g(ica)%MS(i,n)%dcondt(j)*acc_mod(i,j,ica,n)
      end do
      end do
      if(ica/=2) then
        j=3
!        do k = 1, g(ica)%n_masscom
!          do in=1,g(ica)%n_bin*g(ica)%L
!            n=(in-1)/g(ica)%N_BIN+1
!            i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
        do k = 1, g(ica)%n_masscom
           g(ica)%MS(i,n)%dmassdt(1+k,j) = g(ica)%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,ica,n)
        end do
        end do
        end do
!        do in=1,g(ica)%n_bin*g(ica)%L
!          n=(in-1)/g(ica)%N_BIN+1
!          i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
          g(ica)%MS(i,n)%dcondt(j) = g(ica)%MS(i,n)%dcondt(j)*acc_mod(i,j,ica,n)
        end do
        end do
      else
        j=7
!        do k = 1, g(ica)%n_masscom
!          do in=1,g(ica)%n_bin*g(ica)%L
!            n=(in-1)/g(ica)%N_BIN+1
!            i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
        do k = 1, g(ica)%n_masscom
           g(ica)%MS(i,n)%dmassdt(1+k,j) = g(ica)%MS(i,n)%dmassdt(1+k,j)*acc_mod(i,j,ica,n)
        end do
        end do
        end do
!        do in=1,g(ica)%n_bin*g(ica)%L
!          n=(in-1)/g(ica)%N_BIN+1
!          i=in-(n-1)*g(ica)%N_BIN
        do n = 1, g(ica)%L
        do i = 1, g(ica)%N_BIN
          g(ica)%MS(i,n)%dcondt(j) = g(ica)%MS(i,n)%dcondt(j)*acc_mod(i,j,ica,n)
        end do
        end do
      endif
    case default
       LOG_ERROR("mod_other_tendency_ap_vap",*) "Not ready yet",switch
       call PRC_abort
    end select

  end subroutine mod_other_tendency_ap_vap

end module mod_amps_check
