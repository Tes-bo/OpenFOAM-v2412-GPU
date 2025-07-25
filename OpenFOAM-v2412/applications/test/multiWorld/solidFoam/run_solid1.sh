#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#------------------------------------------------------------------------------

solidFoam -case solid1 -world solid1 2>&1 | tee log.solidFoam_solid1
read dummy

#------------------------------------------------------------------------------
