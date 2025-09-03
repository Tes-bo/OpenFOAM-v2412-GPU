divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Description
#     Driver-specific m4/lemon macros for field expressions.
#
#------------------------------------------------------------------------------

include(`m4/lemon/base-setup.m4')dnl
include([m4/lemon/operator-precedence.m4])dnl
dnl
include([m4/lemon/rules-standard.m4])dnl
include([m4/lemon/rules-operations.m4])dnl
include([m4/lemon/rules-functions.m4])dnl
include([m4/lemon/rules-components.m4])dnl
include([m4/lemon/rules-fields-components.m4])dnl
include([m4/lemon/rules-scalar-logic.m4])dnl
dnl
divert(-1)dnl

use_bool_logic()dnl     # Use boolField directly


#-------------------------------------------------------------------------------
# Driver rules
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Standard rules for fields: declaration, new/get, driver functions etc.

include([m4/lemon/rules-fields.m4])dnl
divert(-1)dnl


#------------------------------------------------------------------------------

# Additional safety measures

undefine([substr])dnl   # Avoid collision with C/C++ naming

#------------------------------------------------------------------------------
divert(0)dnl
