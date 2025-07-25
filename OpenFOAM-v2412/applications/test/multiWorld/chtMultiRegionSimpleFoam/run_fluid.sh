#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#------------------------------------------------------------------------------

chtMultiRegionSimpleFoam -case fluid -world fluid 2>&1 | tee run_fluid.log
read dummy

#------------------------------------------------------------------------------
