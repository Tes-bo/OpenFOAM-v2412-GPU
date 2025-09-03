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
#     doxyFilter-top.awk
#
# Description
#     Only output the first /* ... */ comment section found in the file
#     Use @cond / @endcond to suppress documenting all classes/variables
#     - This is useful for application files in which only the first
#       block documents the application itself.
#
#------------------------------------------------------------------------------

BEGIN {
    state = 0
}

# A '/*' at the beginning of a line starts a comment block
/^ *\/\*/ {
   state++
}

# Check first line
# either started with a comment or skip documentation for the whole file
FNR == 1 {
   if (!state)
   {
      print "//! @cond OpenFOAMIgnoreAppDoxygen"
      state = 2
   }
}

# A '*/' ends the comment block
# skip documentation for rest of the file
/\*\// {
    if (state == 1)
    {
        print
        print "//! @cond OpenFOAMIgnoreAppDoxygen"
    }
    state = 2
    next
}

# Print everything within the first comment block
{
    if (state)
    {
        print
    }
    next
}

END {
    if (state == 2)
    {
        print "//! @endcond"
    }
}

#------------------------------------------------------------------------------
