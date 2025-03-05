Module par_amps
  implicit none
  public
  save
  ! parameters
  !---------------------------------------------
  ! ice particles
  ! array specfication for qipv
  ! NOTE:
  !      - aerosol mass components have to be in consecutive order.
  !
  !      - volume, lengths components have to be in consecutive order
  !        following the concentration.
  !
  !    These are redefined at the initialization based on the configuration.
  !
  !    total mass of ice particle
  integer :: imt_q=1
  !    rime, aggregate, crystal, melt water, frozen water (nucleation)
  integer :: imr_q=10,ima_q=15,imc_q=11,imw_q=12,imf_q=16
  !    total aerosol mass, soluble aerosol mass
  integer :: imat_q=13,imas_q=14
  !    number concentration
  integer :: icon_q=2
  !    circiumcsribing volume (dry), a axis length, c axis length, dendritic length
  integer :: ivcs_q=3,iacr_q=4,iccr_q=5,idcr_q=6
  !    center of gravity coordinates: a and c, and extra crystalline structure
  integer :: iag_q=7,icg_q=8,inex_q=9

  !
  ! array specification for g%MS%mass,g%MS%dmassdt,new_M
  integer :: imt=1
  integer :: imr=2,ima=3,imc=4,imw=8,imf=9, imat=5,imas=6,imai=7
  ! array specification for ratio_M
  integer :: imr_m=1,ima_m=2,imc_m=3,imw_m=7,imf_m=8, imat_m=4,imas_m=5,imai_m=6
  ! array specification for new_Q
  integer :: ivcs=1,iacr=2,iccr=3,idcr=4,iag=5,icg=6,inex=7

  !---------------------------------------------
  ! array specfication for qrpv
  integer :: rmt_q=1
  integer :: rmat_q=3,rmas_q=4
  integer :: rcon_q=2
  ! array specification for g%MS%mass,g%MS%dmassdt,new_M
  integer :: rmt=1
  integer :: rmat=2,rmas=3,rmai=4
  ! array specification for ratio_M
  integer :: rmat_m=1,rmas_m=2,rmai_m=3

  !---------------------------------------------
  ! array specfication for qapv
  integer :: amt_q=1
  integer :: acon_q=2
  integer :: ams_q=3
  ! array specification for g%MS%mass,g%MS%dmassdt,new_M
  integer :: amt=1
  integer :: ams=2,ami=3

  ! number of mass component, starding and ending indexes
  integer :: nvar_mcp_ice, i1_mcp_ice, i2_mcp_ice
  ! number of volumne component
  integer :: nvar_vcp_ice, i1_vcp_ice, i2_vcp_ice
  ! number of axis component, starting and ending indexes
  integer :: nvar_acp_ice, i1_acp_ice, i2_acp_ice
  ! number of concentration component, starting and ending indexes
  integer :: nvar_ccp_ice, i1_ccp_ice, i2_ccp_ice
  ! number of non-mass component
  integer :: nvar_nonmcp_ice
  ! number of property variables
  integer :: nvar_pv_ice

  ! number of mass component, starding and ending indexes
  integer :: nvar_mcp_liq, i1_mcp_liq, i2_mcp_liq
  ! number of concentration component, starting and ending indexes
  integer :: nvar_ccp_liq, i1_ccp_liq, i2_ccp_liq
  ! number of property variables
  integer :: nvar_pv_liq
  ! number of non-mass component
  integer :: nvar_nonmcp_liq

  ! number of mass component, starding and ending indexes
  integer :: nvar_mcp_aer, i1_mcp_aer, i2_mcp_aer
  ! number of property variables
  integer :: nvar_pv_aer
  ! number of concentration component, starting and ending indexes
  integer :: nvar_ccp_aer, i1_ccp_aer, i2_ccp_aer


  !---------------------------------------------------
  ! largest bin number that splits cloud and precipitating particles
  integer :: isplit_bin_liq=10
  integer :: isplit_bin_ice=0


  !---------------------------------------------------
  ! debug timing
  integer :: bug_time

  ! debug layer
  integer :: bug_layer

end Module par_amps
