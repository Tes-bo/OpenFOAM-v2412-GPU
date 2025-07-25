#!/bin/bash
. $WM_PROJECT_DIR/etc/bashrc
(cd solid && FOAM_ABORT=true solidFoam -world solid >& log.solidFoam)
