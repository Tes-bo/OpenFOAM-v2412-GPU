#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# File
#     config.csh/example/prefs.csh
#
# Description
#     Example of preset variables for configuring OpenFOAM (C-shell)
#
#     Copy to OpenFOAM-*/etc (or ~/.OpenFOAM) and it will be sourced by
#     OpenFOAM-*/etc/cshrc
#
# See also
#     'foamEtcFile -help' or 'foamEtcFile -list' for the paths searched
#
#------------------------------------------------------------------------------

setenv WM_COMPILER_TYPE ThirdParty
setenv WM_COMPILER Clang
setenv WM_MPLIB SYSTEMOPENMPI

#------------------------------------------------------------------------------
