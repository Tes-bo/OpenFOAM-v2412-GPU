#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#------------------------------------------------------------------------------

chtMultiRegionSimpleFoam -case solid -world solid 2>&1 | tee run_solid.log
read dummy

#------------------------------------------------------------------------------
