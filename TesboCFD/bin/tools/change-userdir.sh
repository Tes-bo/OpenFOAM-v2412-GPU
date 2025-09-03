#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2017-2022 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Script
#     . change-userdir.sh PREFIX [SUFFIX]
#
#     Shortcuts (prefix)
#         -home           "$HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION"
#         -none           remove from environment
#
#     Shortcuts (suffix)
#         -platforms      "platforms/$WM_OPTIONS"
#
# Description
#     Change WM_PROJECT_USER_DIR, FOAM_USER_APPBIN, FOAM_USER_LIBBIN
#     and the respective entries in PATH, LD_LIBRARY_PATH.
#     Also adjusts FOAM_RUN.
#
#     This can be useful with compiling additional OpenFOAM programs
#     (that use FOAM_USER_APPBIN, FOAM_USER_LIBBIN for their build),
#     to avoid conflicts with the normal user bin/lib files.
#
#     The suffix value should normally include "platforms/$WM_OPTIONS"
#
# Example
#     . /path/change-userdir.sh -home -platforms
#
#   corresponds to the standard user location:
#
#     $HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION/platforms/$WM_OPTIONS
#
#------------------------------------------------------------------------------

if [ "$#" -ge 1 ]
then
    prefix="$1"
    suffix="$2"

    foamOldDirs="$FOAM_USER_APPBIN $FOAM_USER_LIBBIN"
    foamClean="$WM_PROJECT_DIR/bin/foamCleanPath"
    if [ -x "$foamClean" ]
    then
        cleaned=$($foamClean -env=PATH "$foamOldDirs") && PATH="$cleaned"
        cleaned=$($foamClean -env=LD_LIBRARY_PATH "$foamOldDirs") \
            && LD_LIBRARY_PATH="$cleaned"
    fi

    case "$suffix" in
        -plat*) suffix="platforms/$WM_OPTIONS" ;;
    esac
    case "$prefix" in
        -home)  prefix="$HOME/OpenFOAM/$USER-${WM_PROJECT_VERSION:-unknown}" ;;
        -none)  unset prefix ;;
    esac

    if [ -n "$prefix" ]
    then
        export WM_PROJECT_USER_DIR="$prefix"
        export FOAM_RUN="$prefix/run"

        prefix="$prefix${suffix:+/}${suffix}"
        export FOAM_USER_APPBIN="$prefix/bin"
        export FOAM_USER_LIBBIN="$prefix/lib"

        PATH="$FOAM_USER_APPBIN:$PATH"
        LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"
    else
        unset WM_PROJECT_USER_DIR FOAM_RUN FOAM_USER_APPBIN FOAM_USER_LIBBIN
    fi
fi

unset foamClean foamOldDirs cleaned prefix suffix

#------------------------------------------------------------------------------
