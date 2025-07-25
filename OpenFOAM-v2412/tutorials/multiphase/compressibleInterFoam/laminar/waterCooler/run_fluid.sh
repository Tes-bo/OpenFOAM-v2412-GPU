#!/bin/bash
. $WM_PROJECT_DIR/etc/bashrc
(cd fluid && FOAM_ABORT=true compressibleInterFoam -world fluid >& log.compressibleInterFoam)
