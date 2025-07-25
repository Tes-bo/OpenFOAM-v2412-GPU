#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2017-2020 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Script
#     . change-sitedir.sh PREFIX [SUFFIX]
#
#     Shortcuts (prefix)
#         -prefix         "$WM_PROJECT_DIR/../site"
#         -project        "$WM_PROJECT_DIR/site"
#         -none           remove from environment
#
#     Shortcuts (suffix)
#         -platforms      "platforms/$WM_OPTIONS"
#
# Description
#     Change WM_PROJECT_SITE, FOAM_SITE_APPBIN, FOAM_SITE_LIBBIN
#     and the respective entries in PATH, LD_LIBRARY_PATH.
#
#     This can be useful when temporarily reassigning the site directory
#     when packaging OpenFOAM.
#
#     The suffix value should normally include "platforms/$WM_OPTIONS"
#
# Example
#     . /path/change-sitedir.sh -project -platforms
#
#   corresponds to the standard site location:
#
#     $WM_PROJECT_DIR/site{/$FOAM_API/platforms/$WM_OPTIONS}
#
#------------------------------------------------------------------------------

if [ "$#" -ge 1 ]
then
    prefix="$1"
    suffix="$2"

    foamOldDirs="$FOAM_SITE_APPBIN $FOAM_SITE_LIBBIN \
        $WM_PROJECT_SITE $WM_PROJECT_DIR/site"
    foamClean="$WM_PROJECT_DIR/bin/foamCleanPath"
    if [ -x "$foamClean" ]
    then
        cleaned=$($foamClean -env=PATH "$foamOldDirs") && PATH="$cleaned"
        cleaned=$($foamClean -env=LD_LIBRARY_PATH "$foamOldDirs") \
            && LD_LIBRARY_PATH="$cleaned"
    fi

    case "$suffix" in
        -plat*)     suffix="platforms/$WM_OPTIONS" ;;
    esac
    case "$prefix" in
        -prefix)    prefix="${WM_PROJECT_DIR%/*}/site" ;;
        -project)   prefix="$WM_PROJECT_DIR/site" ;;
        -none)      unset prefix ;;
    esac

    if [ -n "$prefix" ]
    then
        export WM_PROJECT_SITE="$prefix"

        prefix="$prefix/${WM_PROJECT_VERSION:-unknown}${suffix:+/}${suffix}"

        export FOAM_SITE_APPBIN="$prefix/bin"
        export FOAM_SITE_LIBBIN="$prefix/lib"
        PATH="$FOAM_SITE_APPBIN:$PATH"
        LD_LIBRARY_PATH="$FOAM_SITE_LIBBIN:$LD_LIBRARY_PATH"
    else
        unset WM_PROJECT_SITE FOAM_SITE_APPBIN FOAM_SITE_LIBBIN
    fi
fi

unset foamClean foamOldDirs cleaned prefix suffix

#------------------------------------------------------------------------------
