#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#     Copyright (C) 2020 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# File
#     config.sh/example/prefs.sh
#
# Description
#     Example of preset variables for configuring OpenFOAM (POSIX shell)
#
#     Copy to OpenFOAM-*/etc (or ~/.OpenFOAM) and it will be sourced by
#     OpenFOAM-*/etc/bashrc
#
# See also
#     'foamEtcFile -help' or 'foamEtcFile -list' for the paths searched
#
#------------------------------------------------------------------------------

export WM_COMPILER_TYPE=ThirdParty
export WM_COMPILER=Clang
export WM_MPLIB=SYSTEMOPENMPI

#------------------------------------------------------------------------------
