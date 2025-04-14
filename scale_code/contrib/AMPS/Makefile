################################################################################
#
# Makefile for scale library
#
################################################################################

TOPDIR     = ../..

SYSDEP_DIR = $(TOPDIR)/sysdep
include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

SCALE_BUILD_DIR ?= $(SCALELIBDIR)/src/.libs$(POSTFIX)

BUILD_DIR = .libs$(POSTFIX)
LIBNAME = libamps$(POSTFIX).a

OBJS = \
	$(BUILD_DIR)/acc_amps.o	\
	$(BUILD_DIR)/par_amps.o	\
	$(BUILD_DIR)/com_amps.o	\
	$(BUILD_DIR)/maxdims.o	\
	$(BUILD_DIR)/mod_amps_const.o	\
	$(BUILD_DIR)/mod_scladv_slppm.o	\
	$(BUILD_DIR)/mod_amps_utility.o	\
	$(BUILD_DIR)/class_Thermo_Var.o	\
	$(BUILD_DIR)/class_AirGroup.o	\
	$(BUILD_DIR)/class_Ice_Shape.o	\
	$(BUILD_DIR)/class_Mass_Bin.o	\
	$(BUILD_DIR)/class_Group.o	\
	$(BUILD_DIR)/mod_amps_core.o	\
	$(BUILD_DIR)/mod_amps_check.o	\
	$(BUILD_DIR)/class_Cloud_Micro.o	\
	$(BUILD_DIR)/mod_amps_lib.o


all:
	@echo;echo "Entering AMPS."
	$(MAKE) makedir
	$(MAKE) $(LIBNAME)
	@echo "Complete making AMPS."

makedir:
	mkdir -p $(BUILD_DIR)

$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -f $(OBJS) $(BUILD_DIR)/*.mod *.lst *~

distclean: clean
	rm -f $(LIBNAME)

allclean: distclean
	rm -rf $(BUILD_DIR)
	rm -rf .libs*
	rm -f *.a


.SUFFIXES:
.SUFFIXES: .o .F90

$(BUILD_DIR)/%.o : %.F90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) $(MODDIROPT) $(BUILD_DIR) -I$(SCALELIBDIR)/include -I$(SCALE_BUILD_DIR) -I$(BUILD_DIR) -o $@ -c $<

.PHONY : clean


ifneq ($(SCALE_IGNORE_SRCDEP),T)

#$(BUILD_DIR)/acc_amps.o:
#$(BUILD_DIR)/par_amps.o:
$(BUILD_DIR)/com_amps.o: $(BUILD_DIR)/acc_amps.o
$(BUILD_DIR)/maxdims.o: $(BUILD_DIR)/acc_amps.o
$(BUILD_DIR)/mod_amps_const.o: $(BUILD_DIR)/acc_amps.o
$(BUILD_DIR)/mod_scladv_slppm.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/com_amps.o
$(BUILD_DIR)/mod_amps_utility.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_scladv_slppm.o
$(BUILD_DIR)/class_Thermo_Var.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o
$(BUILD_DIR)/class_AirGroup.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o $(BUILD_DIR)/class_Thermo_Var.o
$(BUILD_DIR)/class_Ice_Shape.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o $(BUILD_DIR)/class_Thermo_Var.o
$(BUILD_DIR)/class_Mass_Bin.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o $(BUILD_DIR)/class_Thermo_Var.o $(BUILD_DIR)/class_Ice_Shape.o
$(BUILD_DIR)/class_Group.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o $(BUILD_DIR)/class_Thermo_Var.o $(BUILD_DIR)/class_AirGroup.o $(BUILD_DIR)/class_Ice_Shape.o $(BUILD_DIR)/class_Mass_Bin.o
$(BUILD_DIR)/mod_amps_core.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/class_Thermo_Var.o $(BUILD_DIR)/class_AirGroup.o $(BUILD_DIR)/class_Group.o
$(BUILD_DIR)/mod_amps_check.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/class_Thermo_Var.o $(BUILD_DIR)/class_AirGroup.o $(BUILD_DIR)/class_Group.o
$(BUILD_DIR)/class_Cloud_Micro.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o $(BUILD_DIR)/class_AirGroup.o $(BUILD_DIR)/class_Mass_Bin.o $(BUILD_DIR)/class_Group.o $(BUILD_DIR)/mod_amps_core.o $(BUILD_DIR)/mod_amps_check.o
$(BUILD_DIR)/mod_amps_lib.o: $(BUILD_DIR)/acc_amps.o $(BUILD_DIR)/par_amps.o $(BUILD_DIR)/com_amps.o $(BUILD_DIR)/maxdims.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/mod_amps_utility.o $(BUILD_DIR)/mod_amps_const.o $(BUILD_DIR)/class_AirGroup.o $(BUILD_DIR)/class_Group.o $(BUILD_DIR)/class_Cloud_Micro.o

endif
