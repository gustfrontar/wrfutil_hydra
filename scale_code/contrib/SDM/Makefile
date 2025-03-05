#################################################################
#
# Makefile for Super-Droplet Method
#
#################################################################

TOPDIR      = ../..

SYSDEP_DIR = $(TOPDIR)/sysdep
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

SCALE_BUILD_DIR ?= $(SCALELIBDIR)/src/.libs$(POSTFIX)

BUILD_DIR = .libs$(POSTFIX)
LIBNAME = libsdm$(POSTFIX).a

DEPENDLIB =

MODS	= \
	  $(BUILD_DIR)/m_sdm_calc.mod \
	  $(BUILD_DIR)/m_sdm_condensation_water.mod \
	  $(BUILD_DIR)/m_sdm_meltfreeze.mod \
	  $(BUILD_DIR)/m_sdm_subldep.mod \
	  $(BUILD_DIR)/m_sdm_sd2fluid.mod \
	  $(BUILD_DIR)/m_sdm_coalescence.mod \
	  $(BUILD_DIR)/m_sdm_coalescence_cold.mod \
	  $(BUILD_DIR)/m_sdm_idutil.mod \
	  $(BUILD_DIR)/m_sdm_motion.mod \
	  $(BUILD_DIR)/m_sdm_boundary.mod \
	  $(BUILD_DIR)/m_sdm_fluidconv.mod \
	  $(BUILD_DIR)/m_sdm_io.mod \
	  $(BUILD_DIR)/m_sdm_coordtrans.mod \
	  $(BUILD_DIR)/m_sdm_memmgr.mod \
	  $(BUILD_DIR)/m_sdm_ptlformation.mod \
	  $(BUILD_DIR)/m_sdm_common.mod

OBJS	= \
	  $(BUILD_DIR)/sdm_calc.o \
	  $(BUILD_DIR)/sdm_condensation_water.o \
	  $(BUILD_DIR)/sdm_meltfreeze.o \
	  $(BUILD_DIR)/sdm_subldep.o \
	  $(BUILD_DIR)/sdm_sd2fluid.o \
	  $(BUILD_DIR)/sdm_coalescence.o \
	  $(BUILD_DIR)/sdm_coalescence_cold.o \
	  $(BUILD_DIR)/sdm_idutil.o \
	  $(BUILD_DIR)/sdm_motion.o \
	  $(BUILD_DIR)/sdm_boundary.o \
	  $(BUILD_DIR)/sdm_fluidconv.o \
	  $(BUILD_DIR)/sdm_io.o \
	  $(BUILD_DIR)/sdm_coordtrans.o \
	  $(BUILD_DIR)/sdm_memmgr.o \
	  $(BUILD_DIR)/sdm_ptlformation.o \
	  $(BUILD_DIR)/sdm_common.o

all:
	@echo;echo "Entering contrib/SDM..."
	$(MAKE) -C $(MTDIR)
	$(MAKE) -C $(GAGETRYDIR)
	$(MAKE) $(LIBNAME)

makedir:
	mkdir -p $(BUILD_DIR)

$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

$(BUILD_DIR)/sdm_calc.o: sdm_calc.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_fluidconv.o \
	$(BUILD_DIR)/sdm_sd2fluid.o \
	$(BUILD_DIR)/sdm_boundary.o \
	$(BUILD_DIR)/sdm_motion.o \
	$(BUILD_DIR)/sdm_coalescence.o \
	$(BUILD_DIR)/sdm_coalescence_cold.o \
	$(BUILD_DIR)/sdm_condensation_water.o \
	$(BUILD_DIR)/sdm_meltfreeze.o \
	$(BUILD_DIR)/sdm_subldep.o \
	$(BUILD_DIR)/sdm_ptlformation.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_prc.o \
	$(SCALE_BUILD_DIR)/scale_time.o \
	$(SCALE_BUILD_DIR)/scale_comm_cartesC.o \
	$(SCALE_BUILD_DIR)/scale_const.o \

$(BUILD_DIR)/sdm_subldep.o: sdm_subldep.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_fluidconv.o \
	$(SCALE_BUILD_DIR)/scale_atmos_hydrometeor.o \
	$(SCALE_BUILD_DIR)/scale_atmos_saturation.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_precision.o

$(BUILD_DIR)/sdm_meltfreeze.o: sdm_meltfreeze.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_fluidconv.o \
	$(SCALE_BUILD_DIR)/scale_atmos_hydrometeor.o \
	$(SCALE_BUILD_DIR)/scale_atmos_saturation.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_precision.o

$(BUILD_DIR)/sdm_condensation_water.o: sdm_condensation_water.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_fluidconv.o \
	$(SCALE_BUILD_DIR)/scale_atmos_hydrometeor.o \
	$(SCALE_BUILD_DIR)/scale_time.o \
	$(SCALE_BUILD_DIR)/scale_atmos_saturation.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_precision.o

$(BUILD_DIR)/sdm_sd2fluid.o: sdm_sd2fluid.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_precision.o

$(BUILD_DIR)/sdm_coalescence.o: sdm_coalescence.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_idutil.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(MTDIR)/rng_uniform_mt.mod \
	$(GAGETRYDIR)/gadg_algorithm.mod

$(BUILD_DIR)/sdm_coalescence_cold.o: sdm_coalescence_cold.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_idutil.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_atmos_saturation.o \
	$(MTDIR)/rng_uniform_mt.mod \
	$(GAGETRYDIR)/gadg_algorithm.mod

$(BUILD_DIR)/sdm_idutil.o: sdm_idutil.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_prc.o \
	$(GAGETRYDIR)/gadg_algorithm.mod

$(BUILD_DIR)/sdm_motion.o: sdm_motion.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC_real.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_prc.o

$(BUILD_DIR)/sdm_boundary.o: sdm_boundary.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_common.o \
	$(SCALE_BUILD_DIR)/scale_prc.o \
	$(SCALE_BUILD_DIR)/scale_io.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC.o

$(BUILD_DIR)/sdm_fluidconv.o: sdm_fluidconv.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_atmos_thermodyn.o \
	$(SCALE_BUILD_DIR)/scale_precision.o

$(BUILD_DIR)/sdm_io.o: sdm_io.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(SCALE_BUILD_DIR)/scale_time.o \
	$(SCALE_BUILD_DIR)/scale_prc.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_io.o

$(BUILD_DIR)/sdm_coordtrans.o: sdm_coordtrans.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(SCALE_BUILD_DIR)/scale_io.o \
	$(SCALE_BUILD_DIR)/scale_prc.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC_real.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC.o \
	$(SCALE_BUILD_DIR)/scale_precision.o

$(BUILD_DIR)/sdm_memmgr.o: sdm_memmgr.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o

$(BUILD_DIR)/sdm_ptlformation.o: sdm_ptlformation.f90 \
	$(DEPENDLIB) \
	$(BUILD_DIR)/sdm_common.o \
	$(BUILD_DIR)/sdm_coordtrans.o \
	$(BUILD_DIR)/sdm_fluidconv.o \
	$(BUILD_DIR)/sdm_sd2fluid.o \
	$(BUILD_DIR)/sdm_condensation_water.o \
	$(BUILD_DIR)/sdm_meltfreeze.o \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(SCALE_BUILD_DIR)/scale_prc.o \
	$(SCALE_BUILD_DIR)/scale_prc_cartesC.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC.o \
	$(SCALE_BUILD_DIR)/scale_io.o \
	$(SCALE_BUILD_DIR)/scale_topography.o \
	$(SCALE_BUILD_DIR)/scale_atmos_grid_cartesC_real.o

$(BUILD_DIR)/sdm_common.o: sdm_common.f90 \
	$(DEPENDLIB) \
	$(SCALE_BUILD_DIR)/scale_precision.o \
	$(SCALE_BUILD_DIR)/scale_io.o \
	$(SCALE_BUILD_DIR)/scale_const.o \
	$(MTDIR)/rng_uniform_mt.mod

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.lst
	$(MAKE) -C $(GAGETRYDIR) clean
	$(MAKE) -C $(MTDIR) clean

distclean: clean
	rm -f $(LIBNAME)
	$(MAKE) -C $(GAGETRYDIR) distclean
	$(MAKE) -C $(MTDIR) distclean

allclean: distclean
	rm -f *.a
	rm -rf .libs*
	$(MAKE) -C $(GAGETRYDIR) allclean
	$(MAKE) -C $(MTDIR) allclean

.SUFFIXES:
.SUFFIXES: .o .f90 .c

$(BUILD_DIR)/%.o: %.f90
	$(FC) $(FFLAGS) $(FFLAGS_SDM) $(ADDITIONAL_FFLAGS) $(MODDIROPT) $(BUILD_DIR) -I$(SCALELIBDIR)/include -I$(SCALE_BUILD_DIR) -I$(BUILD_DIR) -I$(MTDIR) -I$(GAGETRYDIR) $(SCALE_NETCDF_INCLUDE) -o $@ -c $<

%.o: %.mod

.PHONY : clean distclean allclean
