#------------------------------------------------------------------------------
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
# Script
#     doxyFilter-ignore.awk
#
# Description
#     - Prefix file contents with doxygen @file tag and %filePath% tag
#       that will be changed in a subsequent sed script
#     - Surround the contents of an entire file with @cond / @endcond
#       to skip documenting all classes/variables
#
#------------------------------------------------------------------------------

BEGIN {
   print "//! @file %filePath%"
   print "//! @cond OpenFOAMIgnoreAppDoxygen"
}

{ print }

END {
   print "//! @endcond"
}

#------------------------------------------------------------------------------
