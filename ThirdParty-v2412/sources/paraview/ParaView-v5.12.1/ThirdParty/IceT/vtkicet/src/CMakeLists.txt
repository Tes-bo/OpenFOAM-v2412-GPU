## Copyright 2003 Sandia Coporation
## Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
## the U.S. Government retains certain rights in this software.
##
## This source code is released under the New BSD License.
#

ADD_SUBDIRECTORY(ice-t)
ADD_SUBDIRECTORY(communication)
if (FALSE) # XXX(kitware): ParaView doesn't need this.
ADD_SUBDIRECTORY(strategies)
endif ()

IF (ICET_USE_OPENGL)
  ADD_SUBDIRECTORY(gl)
ENDIF (ICET_USE_OPENGL)
