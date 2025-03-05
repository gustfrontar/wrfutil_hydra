################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR     = $(abspath ../..)
SYSDEP_DIR = $(abspath ../../sysdep)

BUILD_DIR = "undefined"
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

LIBNAME = libscale$(POSTFIX).a

DEPENDLIB =

VPATH = common:                          \
		common/fem:                      \
        io:                              \
        file:                            \
        communication:                   \
        atmosphere/common:               \
        atmosphere/grid:                 \
        atmosphere/physics/aerosol:      \
        atmosphere/physics/chemistry:    \
        atmosphere/physics/cumulus:      \
        atmosphere/physics/lightning:    \
        atmosphere/physics/microphysics: \
        atmosphere/physics/pbl:          \
        atmosphere/physics/radiation:    \
        atmosphere/physics/surface:      \
        land:                            \
        land/grid:                       \
        land/dynamics:                   \
        land/physics:                    \
        ocean:                           \
        ocean/grid:                      \
        ocean/dynamics:                  \
        ocean/physics:                   \
        urban:                           \
        urban/grid:                      \
        urban/dynamics:                  \
        coupler:                         \
        coupler/physics/surface:         \
        data-assim:                      \
        \
        atmos-rm/dynamics:               \
        atmos-rm/turbulence:

LIBVERSION = $(shell git rev-parse --short HEAD 2> /dev/null)
ifeq ($(LIBVERSION),)
  LIBVERSION  = $(shell cat VERSION)
else
  LIBVERSION := $(LIBVERSION)
endif

SDM_EXIST := $(shell find $(TOPDIR)/contrib -name SDM 2> /dev/null )

OBJS_CORE = \
	$(BUILD_DIR)/scale_precision.o                  \
	$(BUILD_DIR)/scale_const.o                      \
	$(BUILD_DIR)/scale_io.o                         \
	$(BUILD_DIR)/scale_time.o                       \
	$(BUILD_DIR)/scale_topography.o                 \
	$(BUILD_DIR)/scale_prc_cartesC.o                \
	$(BUILD_DIR)/scale_comm_cartesC.o               \
	$(BUILD_DIR)/scale_atmos_thermodyn.o            \
	$(BUILD_DIR)/scale_atmos_saturation.o           \
	$(BUILD_DIR)/scale_atmos_grid_cartesC_index.o   \
	$(BUILD_DIR)/scale_atmos_grid_cartesC.o         \
	$(BUILD_DIR)/scale_atmos_grid_cartesC_real.o    \
	$(BUILD_DIR)/scale_atmos_grid_cartesC_metric.o

OBJS = \
	$(BUILD_DIR)/scale.o	\
	\
	$(BUILD_DIR)/scale_sys_sigvars.o \
	$(BUILD_DIR)/scale_file_netcdf.o \
	$(BUILD_DIR)/scale_file_grads.o	\
	\
	$(BUILD_DIR)/scale_calendar.o	\
	$(BUILD_DIR)/scale_debug.o \
	$(BUILD_DIR)/scale_fpm.o \
	$(BUILD_DIR)/scale_index.o \
	$(BUILD_DIR)/scale_interp.o	\
	$(BUILD_DIR)/scale_interp_vert.o \
	$(BUILD_DIR)/scale_mapprojection.o \
	$(BUILD_DIR)/scale_coriolis.o	\
	$(BUILD_DIR)/scale_matrix.o	\
	$(BUILD_DIR)/scale_gmres.o	\
	$(BUILD_DIR)/scale_sort.o	\
	$(BUILD_DIR)/scale_hash.o	\
	$(BUILD_DIR)/scale_filter.o	\
	$(BUILD_DIR)/scale_prc.o	\
	$(BUILD_DIR)/scale_prc_icoA.o	\
	$(BUILD_DIR)/scale_prof.o	\
	$(BUILD_DIR)/scale_random.o	\
	$(BUILD_DIR)/scale_sigvars.o \
	$(BUILD_DIR)/scale_specfunc.o	\
	$(BUILD_DIR)/scale_statistics.o	\
	$(BUILD_DIR)/scale_tracer.o	\
	$(BUILD_DIR)/scale_vector.o	\
	$(BUILD_DIR)/scale_io.o	\
	$(BUILD_DIR)/scale_monitor.o	\
	$(BUILD_DIR)/scale_monitor_cartesC.o \
	$(BUILD_DIR)/scale_comm_cartesC_nest.o \
	$(BUILD_DIR)/scale_comm_ensemble.o \
	$(BUILD_DIR)/scale_comm_icoA.o	\
	$(BUILD_DIR)/scale_file.o	\
	$(BUILD_DIR)/scale_file_cartesC.o	\
	$(BUILD_DIR)/scale_file_external_input.o \
	$(BUILD_DIR)/scale_file_external_input_cartesC.o \
	$(BUILD_DIR)/scale_file_h.o	\
	$(BUILD_DIR)/scale_file_history.o	\
	$(BUILD_DIR)/scale_file_history_cartesC.o	\
	$(BUILD_DIR)/scale_file_tiledata.o \
	$(BUILD_DIR)/scale_timeint_rk_butcher_tab.o \
	$(BUILD_DIR)/scale_timeint_rk.o \
	$(BUILD_DIR)/scale_polynominal.o \
	$(BUILD_DIR)/scale_sparsemat.o \
	$(BUILD_DIR)/scale_fem_element_base.o \
	$(BUILD_DIR)/scale_fem_element_line.o \
	$(BUILD_DIR)/scale_fem_element_quadrilateral.o \
	$(BUILD_DIR)/scale_fem_element_hexahedral.o \
	$(BUILD_DIR)/scale_fem_element_modalfilter.o \
	$(BUILD_DIR)/scale_fem_localmesh_base.o \
	$(BUILD_DIR)/scale_fem_localmesh_1d.o \
	$(BUILD_DIR)/scale_fem_localmesh_2d.o \
	$(BUILD_DIR)/scale_fem_localmesh_3d.o \
	$(BUILD_DIR)/scale_fem_meshutil_1d.o \
	$(BUILD_DIR)/scale_fem_meshutil_2d.o \
	$(BUILD_DIR)/scale_fem_meshutil_3d.o \
	$(BUILD_DIR)/scale_fem_mesh_bndinfo.o \
	$(BUILD_DIR)/scale_fem_mesh_base.o \
	$(BUILD_DIR)/scale_fem_mesh_base1d.o \
	$(BUILD_DIR)/scale_fem_mesh_base2d.o \
	$(BUILD_DIR)/scale_fem_mesh_base3d.o \
	$(BUILD_DIR)/scale_fem_mesh_rectdom2d.o \
	$(BUILD_DIR)/scale_fem_mesh_cubedom3d.o \
	$(BUILD_DIR)/scale_fem_variableinfo.o \
	$(BUILD_DIR)/scale_fem_localmeshfield_base.o \
	$(BUILD_DIR)/scale_fem_meshfield_base.o \
	$(BUILD_DIR)/scale_comm_fem_meshfield_base.o \
	$(BUILD_DIR)/scale_comm_fem_meshfield_rectdom2d.o \
	$(BUILD_DIR)/scale_comm_fem_meshfield_cubedom3d.o \
	$(BUILD_DIR)/scale_file_fem_common.o \
	$(BUILD_DIR)/scale_file_fem_base.o \
	$(BUILD_DIR)/scale_file_fem_history.o \
	$(BUILD_DIR)/scale_file_fem_restart.o \
	\
	$(BUILD_DIR)/scale_atmos_adiabat.o \
	$(BUILD_DIR)/scale_atmos_aerosol.o \
	$(BUILD_DIR)/scale_atmos_bottom.o	\
	$(BUILD_DIR)/scale_atmos_diagnostic.o	\
	$(BUILD_DIR)/scale_atmos_diagnostic_cartesC.o	\
	$(BUILD_DIR)/scale_atmos_hydrometeor.o \
	$(BUILD_DIR)/scale_atmos_hydrostatic.o \
	$(BUILD_DIR)/scale_atmos_profile.o \
	$(BUILD_DIR)/scale_atmos_refstate.o	\
	$(BUILD_DIR)/scale_atmos_solarins.o	\
	$(BUILD_DIR)/scale_atmos_grid_icoA.o \
	$(BUILD_DIR)/scale_atmos_grid_icoA_index.o \
	$(BUILD_DIR)/scale_atmos_dyn.o \
	$(BUILD_DIR)/scale_atmos_dyn_common.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd2.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd4.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd6.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd8.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud1.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud3.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud3Koren1993.o \
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud5.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud7.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_numfilter.o \
	$(BUILD_DIR)/scale_atmos_dyn_fvm_fct.o \
	$(BUILD_DIR)/scale_atmos_dyn_dgm_vars.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_bnd.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_operator.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_hydrostatic.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_rkcommon.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_large.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_large_euler.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_large_rk3.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk3.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk4.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk7s6o.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk11s8o.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rkdg.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_euler.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_rk3.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_linrk.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_rkdg.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_large.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_large_fvm_heve.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_large_dgm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_fvm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_dgm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_fvm_hevi.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_dgm_hevi.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_dgm_hevi_common.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_fvm_hivi.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_dgm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve_ref.o \
	$(BUILD_DIR)/scale_atmos_phy_ae_kajino13.o \
	$(BUILD_DIR)/scale_atmos_phy_ae_offline.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_admin.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_common.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_aerocarb.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_aerosulf.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_aerosalt.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_aerodust.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_aerotrac.o	\
	$(BUILD_DIR)/scale_atmos_phy_ae_sprintars_ccn.o	\
	$(BUILD_DIR)/scale_atmos_phy_bl_common.o \
	$(BUILD_DIR)/scale_atmos_phy_bl_mynn.o \
	$(BUILD_DIR)/scale_atmos_phy_bl_mynn_jmapplib.o	\
	$(BUILD_DIR)/scale_atmos_phy_ch_rn222.o	\
	$(BUILD_DIR)/scale_atmos_sfc_ch_rn222.o	\
	$(BUILD_DIR)/scale_atmos_phy_cp_common.o \
	$(BUILD_DIR)/scale_atmos_phy_cp_kf.o \
	$(BUILD_DIR)/scale_atmos_phy_cp_kf_jmapplib.o	\
	$(BUILD_DIR)/scale_atmos_phy_lt_sato2019.o \
	$(BUILD_DIR)/scale_atmos_phy_mp_common.o \
	$(BUILD_DIR)/scale_atmos_phy_mp_kessler.o	\
	$(BUILD_DIR)/scale_atmos_phy_mp_sn14.o	\
	$(BUILD_DIR)/scale_atmos_phy_mp_suzuki10.o \
	$(BUILD_DIR)/scale_atmos_phy_mp_tomita08.o \
	$(BUILD_DIR)/scale_atmos_phy_mp_amps.o	\
	$(BUILD_DIR)/scale_atmos_phy_mp_sdm.o \
	$(BUILD_DIR)/scale_atmos_phy_rd_common.o \
	$(BUILD_DIR)/scale_atmos_phy_rd_mstrnx.o \
	$(BUILD_DIR)/scale_atmos_phy_rd_offline.o	\
	$(BUILD_DIR)/scale_atmos_phy_rd_profile.o	\
	$(BUILD_DIR)/scale_atmos_phy_sf_bulk.o \
	$(BUILD_DIR)/scale_atmos_phy_sf_const.o	\
	\
	$(BUILD_DIR)/scale_atmos_phy_tb_common.o \
	$(BUILD_DIR)/scale_atmos_phy_tb_common_dgm.o \
	$(BUILD_DIR)/scale_atmos_phy_tb_dgm_vars.o \
	$(BUILD_DIR)/scale_atmos_phy_tb_d1980.o	\
	$(BUILD_DIR)/scale_atmos_phy_tb_dns.o	\
	$(BUILD_DIR)/scale_atmos_phy_tb_smg.o	\
	$(BUILD_DIR)/scale_atmos_phy_tb_smg_dgm.o	\
	\
	$(BUILD_DIR)/scale_ocean_grid_cartesC.o	\
	$(BUILD_DIR)/scale_ocean_grid_cartesC_index.o	\
	$(BUILD_DIR)/scale_ocean_grid_cartesC_real.o \
	$(BUILD_DIR)/scale_ocean_grid_icoA.o \
	$(BUILD_DIR)/scale_ocean_grid_icoA_index.o \
	$(BUILD_DIR)/scale_ocean_dyn_slab.o	\
	$(BUILD_DIR)/scale_ocean_dyn_offline.o \
	$(BUILD_DIR)/scale_ocean_phy_ice_simple.o	\
	$(BUILD_DIR)/scale_ocean_phy_albedo.o	\
	$(BUILD_DIR)/scale_ocean_phy_albedo_nakajima00.o \
	$(BUILD_DIR)/scale_ocean_phy_roughness.o \
	$(BUILD_DIR)/scale_ocean_phy_roughness_miller92.o	\
	$(BUILD_DIR)/scale_ocean_phy_roughness_moon07.o	\
	$(BUILD_DIR)/scale_ocean_phy_tc.o	\
	$(BUILD_DIR)/scale_land_grid_cartesC.o \
	$(BUILD_DIR)/scale_land_grid_cartesC_index.o \
	$(BUILD_DIR)/scale_land_grid_cartesC_real.o	\
	$(BUILD_DIR)/scale_land_grid_icoA.o	\
	$(BUILD_DIR)/scale_land_grid_icoA_index.o	\
	$(BUILD_DIR)/scale_land_dyn_bucket.o \
	$(BUILD_DIR)/scale_land_phy_matsiro.o	\
	$(BUILD_DIR)/scale_land_phy_snow_diagnos.o \
	$(BUILD_DIR)/scale_land_phy_snow_ky90.o	\
	$(BUILD_DIR)/scale_urban_grid_cartesC.o	\
	$(BUILD_DIR)/scale_urban_grid_cartesC_index.o	\
	$(BUILD_DIR)/scale_urban_grid_cartesC_real.o \
	$(BUILD_DIR)/scale_urban_grid_icoA.o \
	$(BUILD_DIR)/scale_urban_grid_icoA_index.o \
	$(BUILD_DIR)/scale_urban_dyn_kusaka01.o	\
	$(BUILD_DIR)/scale_bulkflux.o	\
	$(BUILD_DIR)/scale_landuse.o \
	$(BUILD_DIR)/scale_cpl_sfc_index.o \
	$(BUILD_DIR)/scale_cpl_phy_sfc_fixed_temp.o	\
	$(BUILD_DIR)/scale_cpl_phy_sfc_skin.o \
	$(BUILD_DIR)/scale_spnudge.o \
	$(BUILD_DIR)/scale_dft.o \
	\
	$(BUILD_DIR)/read_toshiba.o \
	$(BUILD_DIR)/read_toshiba_mpr.o \
	$(BUILD_DIR)/scale_da_read_pawr_toshiba.o \
	$(BUILD_DIR)/scale_da_read_mp_pawr_toshiba.o \
	$(BUILD_DIR)/scale_letkf.o

OBJS_FOR_DYN = \
	$(BUILD_DIR)/scale_atmos_dyn.o \
	$(BUILD_DIR)/scale_atmos_dyn_common.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd2.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd4.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd6.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_cd8.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud1.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud3.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud3Koren1993.o \
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud5.o	\
	$(BUILD_DIR)/scale_atmos_dyn_fvm_flux_ud7.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_vars.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_bnd.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_operator.o	\
	$(BUILD_DIR)/scale_atmos_dyn_dgm_hydrostatic.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_large.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_large_euler.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_large_rk3.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk3.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk4.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rk7s6o.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_short_rkdg.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_euler.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_rk3.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_linrk.o \
	$(BUILD_DIR)/scale_atmos_dyn_tinteg_tracer_rkdg.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_large.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_large_fvm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_large_dgm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short.o \
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_fvm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_dgm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_fvm_hevi.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_short_fvm_hivi.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_dgm_heve.o	\
	$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve_ref.o

MODS = $(addprefix $(BUILD_DIR)/, \
	$(shell echo $(patsubst $(BUILD_DIR)/%.o,%.mod, \
	$(filter-out $(BUILD_DIR)/scale_sys_sigvars.o, \
	$(filter-out $(BUILD_DIR)/scale_file_netcdf.o, \
	$(filter-out $(BUILD_DIR)/read_toshiba.o, \
	$(filter-out $(BUILD_DIR)/read_toshiba_mpr.o, \
	$(OBJS) $(OBJS_CORE)))))) | tr A-Z a-z) )

all:
	$(MAKE) build
	$(MAKE) install

build:
	@echo;echo "Entering scalelib..."; echo "Current version is " $(LIBVERSION)
	$(MAKE) makedir
	$(MAKE) makeobjs
	@echo "Complete making scalelib."

install:
	mkdir -p $(MODDIR)
	mkdir -p $(LIBDIR)
	$(MAKE) modules
	$(MAKE) $(LIBDIR)/$(LIBNAME)

info:
	@$(MAKE) --no-print-directory envlog
	@echo
	@echo "SCALElib"
	@$(MAKE) --no-print-directory conflog


$(LIBDIR)/$(LIBNAME): $(BUILD_DIR)/$(LIBNAME)
	install $< $@

$(BUILD_DIR)/$(LIBNAME): $(OBJS_CORE) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

modules: $(MODS)
	install $? $(MODDIR)

makedir:
	$(MAKE) -C $(AMPSDIR) makedir
	$(MAKE) -C $(SDMDIR) makedir
	mkdir -p $(BUILD_DIR)

ifeq ($(SCALE_IGNORE_SRCDEP),T)
makeobjs: makeamps makesdm $(OBJS)
else
makeobjs: $(OBJS)
endif

$(AMPSDIR)/libamps$(POSTFIX).a: makeamps
makeamps: $(OBJS_CORE)
	$(MAKE) -C $(AMPSDIR) SCALE_BUILD_DIR=$(BUILD_DIR)

$(SDMDIR)/libsdm$(POSTFIX).a: makesdm
makesdm: $(OBJS_CORE)
	$(MAKE) -C $(SDMDIR) SCALE_BUILD_DIR=$(BUILD_DIR)

allclean: distclean
	$(MAKE) -C $(AMPSDIR) allclean
	$(MAKE) -C $(SDMDIR) allclean
	rm -rf $(MODDIR)
	rm -rf $(LIBDIR)

distclean: clean
	rm -f $(LIBDIR)/$(LIBNAME)
	rm -f $(MODDIR)/*.mod
	rm -rf ./.libs*

clean:
	$(MAKE) -C $(AMPSDIR) clean
	$(MAKE) -C $(SDMDIR) clean
	rm -rf $(BUILD_DIR)
	rm -f *.o *.mod *.lst *~

dependency:
	./makedepend .


.SUFFIXES:
.SUFFIXES: .o .F90 .c .erb

%.F90 : %.F90.erb
	erb $< > $@

file/scale_file.F90 : scale_file.F90.erb

file/scale_file_history.F90 : scale_file_history.F90.erb

common/fem/scale_timeint_rk.F90 : scale_timeint_rk.F90.erb

atmos-rm/dynamics/scale_atmos_dyn_fvm_flux.F90 : scale_atmos_dyn_fvm_flux.F90.erb
	erb $< > $@

atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_ud1.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_cd2.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_ud3.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_ud3Koren1993.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_cd4.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_ud5.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_cd6.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_ud7.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@
atmos-rm/dynamics/scale_atmos_dyn_fvm_flux_cd8.F90 : scale_atmos_dyn_fvm_flux_udcd.F90.erb
	fname=$@ erb -T - $< > $@

$(BUILD_DIR)/scale_atmos_dyn_tstep_large_fvm_heve.o: scale_atmos_dyn_tstep_large_fvm_heve.F90
$(BUILD_DIR)/scale_atmos_dyn_tstep_tracer_fvm_heve.o: scale_atmos_dyn_tstep_tracer_fvm_heve.F90

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS)     $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(SCALELIBDIR)/include -I$(BUILD_DIR) $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(CONTRIB_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<
$(OBJS_FOR_DYN) : $(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS_DYN) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(LIBVERSION)\" -I$(SCALELIBDIR)/include -I$(BUILD_DIR) $(SCALE_PAPI_INCLUDE) $(SCALE_NETCDF_INCLUDE) $(CONTRIB_INCLUDE) $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<
$(BUILD_DIR)/%.o : %.c
	$(CC) $(CFLAGS) -I$(SCALELIBDIR)/include $(SCALE_NETCDF_INCLUDE) -o $@ -c $<

.PHONY : clean distclean allclean depend modules dependency

include $(TOPDIR)/utils/make/Make_environments

ifneq ($(SCALE_IGNORE_SRCDEP),T)
-include depend

$(BUILD_DIR)/scale_atmos_phy_mp_amps.o: $(AMPSDIR)/libamps$(POSTFIX).a
$(BUILD_DIR)/scale_atmos_phy_mp_sdm.o: $(SDMDIR)/libsdm$(POSTFIX).a
endif


%.o: %.mod
