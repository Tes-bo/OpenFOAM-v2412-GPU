#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

# 1) Run meshing
# 2) Test input zones and movement

# Meshing
(cd steady && ./Allrun.pre)

if notTest
then
    if canCompile
    then
        (cd code && wmake)
    else
        exit 0
    fi

    . files/RunFunctions

    caseName="movement"

    # Copy/link the steady-state case to movement
    linkParallelCase steady "$caseName"

    # Copy/link support files
    linkFiles files "$caseName"

    # Run
    "$caseName/Allrun.$caseName" $*
fi

#------------------------------------------------------------------------------
