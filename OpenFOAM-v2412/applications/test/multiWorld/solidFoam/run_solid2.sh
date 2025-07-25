#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#------------------------------------------------------------------------------

solidFoam -case solid2 -world solid2 2>&1 | tee log.solidFoam_solid2
read dummy

#------------------------------------------------------------------------------
