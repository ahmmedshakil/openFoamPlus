#----------------------------*- makefile-gmake -*------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Released 2004-2011 OpenCFD Ltd.
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#     Modified code Copyright (C) 2017-2019 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     Makefile
#
# Description
#     Makefile used to compile wmake utilities.
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

# Locally set optimized compilation
WM_COMPILE_OPTION = Opt

GENERAL_RULES   = $(WM_DIR)/rules/General
include $(GENERAL_RULES)/general

archHost    := $(WM_ARCH)$(WM_COMPILER)
archTarget  := $(shell basename $(WMAKE_BIN))


#------------------------------------------------------------------------------
# Targets
#------------------------------------------------------------------------------

.PHONY: all clean message

all: $(WMAKE_BIN)/wmkdepend$(EXT_EXE) message

message:
ifneq ($(archHost),$(archTarget))
	@echo "built wmake-bin ($(archTarget)) for $(archHost) host"
else
	@echo "built wmake-bin ($(archTarget))"
endif

clean:
	@echo "clean wmake-bin ($(archTarget))"
	@rm -rf $(WMAKE_BIN) 2>/dev/null
	@rmdir $(shell dirname $(WMAKE_BIN)) 2>/dev/null || true

$(WMAKE_BIN)/wmkdep$(EXT_EXE): wmkdep.l
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,flex,$(<F))
	$E flex -o $@.c $(<F) && $(cc) $(cFLAGS) $@.c -o $@
	@rm -f $@.c 2>/dev/null

$(WMAKE_BIN)/wmkdepend$(EXT_EXE): wmkdepend.cpp
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,wmkdepend,$(<F))
	$E $(CC) $(c++FLAGS) $(c++LESSWARN) $(<F) -o $@

# $(WMAKE_BIN)/wmkdepend$(EXT_EXE): wmkdepend.rl
# 	@mkdir -p $(WMAKE_BIN)
# 	$(call QUIET_MESSAGE,ragel,$(<F))
# 	$E ragel -G2 -o $@.cpp $(<F) && $(CC) $(c++FLAGS) $(c++LESSWARN) $@.cpp -o $@
# 	@rm -f $@.cpp 2>/dev/null

#------------------------------------------------------------------------------
