#----------------------------*- makefile-gmake -*------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#     Copyright (C) 2017-2020 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# File
#     Makefile
#
# Description
#     Makefile for compiling wmake toolchain.
#
# Note
#     The 'createCode' script is used to pre-generate the scanner code
#     since ragel is not very prevalent.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Use POSIX shell
#------------------------------------------------------------------------------

SHELL   = /bin/sh

#------------------------------------------------------------------------------
# No default suffix rules used
#------------------------------------------------------------------------------

.SUFFIXES:


#------------------------------------------------------------------------------
# Set compilation and dependency building rules
#------------------------------------------------------------------------------

# Locally defined: optimized compilation
WM_COMPILE_OPTION = Opt

# Never want openmp for tools!
WM_COMPILE_CONTROL += ~openmp

GENERAL_RULES   = $(WM_DIR)/rules/General
include $(GENERAL_RULES)/general

archBuild   := $(WM_ARCH)$(WM_COMPILER)
archTarget  := $(shell basename $(WMAKE_BIN))


#------------------------------------------------------------------------------
# Targets
#------------------------------------------------------------------------------

.PHONY: all clean message old

all: $(WMAKE_BIN)/lemon$(EXT_EXE) $(WMAKE_BIN)/wmkdepend$(EXT_EXE) message

# Flex-based processing
old: $(WMAKE_BIN)/wmkdep$(EXT_EXE)

message:
ifneq ("$(archBuild)","$(archTarget)")
	@echo "built wmake-bin ($(archTarget)) with $(archBuild)"
else
	@echo "built wmake-bin ($(archTarget))"
endif

clean:
	@echo "clean wmake-bin ($(archTarget))"
	@rm -rf $(WMAKE_BIN)
	@rmdir $(shell dirname $(WMAKE_BIN)) 2>/dev/null || true

$(WMAKE_BIN)/lemon$(EXT_EXE): lemon.c
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,lemon,$(<F))
	$(cc) $(cFLAGS) $(<F) -o $@

$(WMAKE_BIN)/wmkdep$(EXT_EXE): wmkdep.l
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,flex,$(<F))
	$E flex -o $@.c $(<F) && $(cc) $(cFLAGS) $@.c -o $@
	@rm -f $@.c

$(WMAKE_BIN)/wmkdepend$(EXT_EXE): wmkdepend.cc
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,wmkdepend,$(<F))
	$E $(CC) $(c++FLAGS) $(<F) -o $@

#$(WMAKE_BIN)/wmkdepend$(EXT_EXE): wmkdepend.rl
#	@mkdir -p $(WMAKE_BIN)
#	$(call QUIET_MESSAGE,ragel,$(<F))
#	$E ragel -G2 -o $@.cc $(<F) && $(CC) $(c++FLAGS) $@.cc -o $@
#	@rm -f $@.cc

#------------------------------------------------------------------------------
