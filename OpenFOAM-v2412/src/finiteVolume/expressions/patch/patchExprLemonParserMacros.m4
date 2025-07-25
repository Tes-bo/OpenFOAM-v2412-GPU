divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019-2023 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Description
#     Driver-specific m4/lemon macros for patch expressions.
#
#------------------------------------------------------------------------------

include(`m4/lemon/base-setup.m4')dnl
include([m4/lemon/operator-precedence.m4])dnl
dnl
include([m4/lemon/rules-standard.m4])dnl
include([m4/lemon/rules-operations.m4])dnl
include([m4/lemon/rules-functions.m4])dnl
include([m4/lemon/rules-components.m4])dnl
include([m4/lemon/rules-fields-components.m4])dnl
include([m4/lemon/rules-scalar-logic.m4])dnl
dnl
divert(-1)dnl

use_bool_logic()dnl     # Use boolField directly


#------------------------------------------------------------------------------
# Patch-specific driver rules

#------------------------------------------------------------------------------
# rule_get_patchfields(out, valType, ident)
#
# Description
#     Production rules to get patch internal field
#     and patch neighbour field
#
# Example
#     rule_get_patchfields(sfield, Foam::scalar, SCALAR_ID)
#------------------------------------------------------------------------------

define([rule_get_patchfields],
[dnl
rule_driver_unary_named($1, INTERNAL_FIELD, $3, patchInternalField, $2)dnl
rule_driver_unary_named($1, NEIGHBOUR_FIELD, $3, patchNeighbourField, $2)dnl
rule_driver_unary_named($1, SN_GRAD, $3, patchNormalField, $2)dnl
])


#-------------------------------------------------------------------------------
# Driver rules
#-------------------------------------------------------------------------------

define([rules_driver_surface_functions],
[dnl
_logic_ (lhs) ::= CELL_SET LPAREN identifier (name) RPAREN .dnl
{dnl
    lhs = driver->field_cellSet(make_obj(name)).ptr();
}dnl
_logic_ (lhs) ::= CELL_ZONE LPAREN identifier (name) RPAREN .dnl
{dnl
    lhs = driver->field_cellZone(make_obj(name)).ptr();
}dnl
_logic_ (lhs) ::= FACE_SET LPAREN identifier (name) RPAREN .dnl
{dnl
    lhs = driver->field_faceSet(make_obj(name)).ptr();
}dnl
_logic_ (lhs) ::= FACE_ZONE LPAREN identifier (name) RPAREN .dnl
{dnl
    lhs = driver->field_faceZone(make_obj(name)).ptr();
}dnl
dnl
rule_driver_nullary(_scalar_, FACE_AREA, field_faceArea)dnl
rule_driver_nullary(_vector_, POS, field_faceCentre)dnl  FACE_CENTRE
rule_driver_nullary(_vector_, FACE_EXPR, field_areaNormal)dnl
rule_driver_nullary(_vector_, FACE_NORMAL, field_unitNormal)dnl
dnl
rule_driver_inplace_unary(_scalar_, WEIGHT_AVERAGE, areaAverage)dnl
rule_driver_inplace_unary(_vector_, WEIGHT_AVERAGE, areaAverage)dnl
rule_driver_inplace_unary(_sphTensor_, WEIGHT_AVERAGE, areaAverage)dnl
rule_driver_inplace_unary(_symTensor_, WEIGHT_AVERAGE, areaAverage)dnl
rule_driver_inplace_unary(_tensor_, WEIGHT_AVERAGE, areaAverage)dnl
dnl
rule_driver_inplace_unary(_scalar_, WEIGHT_SUM, areaSum)dnl
rule_driver_inplace_unary(_vector_, WEIGHT_SUM, areaSum)dnl
rule_driver_inplace_unary(_sphTensor_, WEIGHT_SUM, areaSum)dnl
rule_driver_inplace_unary(_symTensor_, WEIGHT_SUM, areaSum)dnl
rule_driver_inplace_unary(_tensor_, WEIGHT_SUM, areaSum)dnl
dnl
])

define([rules_driver_point_functions],
[dnl
rule_driver_nullary(_vector_, POINTS, field_pointField)dnl
dnl
dnl NB use non-driver versions for points - ie, unweighted
dnl
rule_inplace_unary(_scalar_, WEIGHT_AVERAGE, Foam::gAverage)dnl
rule_inplace_unary(_vector_, WEIGHT_AVERAGE, Foam::gAverage)dnl
rule_inplace_unary(_sphTensor_, WEIGHT_AVERAGE, Foam::gAverage)dnl
rule_inplace_unary(_symTensor_, WEIGHT_AVERAGE, Foam::gAverage)dnl
rule_inplace_unary(_tensor_, WEIGHT_AVERAGE, Foam::gAverage)dnl
dnl
rule_inplace_unary(_scalar_, WEIGHT_SUM, Foam::gSum)dnl
rule_inplace_unary(_vector_, WEIGHT_SUM, Foam::gSum)dnl
rule_inplace_unary(_sphTensor_, WEIGHT_SUM, Foam::gSum)dnl
rule_inplace_unary(_symTensor_, WEIGHT_SUM, Foam::gSum)dnl
rule_inplace_unary(_tensor_, WEIGHT_SUM, Foam::gSum)dnl
dnl
])


#------------------------------------------------------------------------------
# rule_faceToPoint(out, in)
# rule_pointToFace(out, in)
#
# Description
#     Production rules for driver faceToPoint, pointToFace,
#     methods
#------------------------------------------------------------------------------

define([rule_faceToPoint],
[rule_driver_unary($1, $2, FACE_TO_POINT, faceToPoint)])

define([rule_pointToFace],
[rule_driver_unary($1, $2, POINT_TO_FACE, pointToFace)])


#------------------------------------------------------------------------------
# Standard rules for fields: declaration, new/get, driver functions etc.

include([m4/lemon/rules-fields.m4])dnl
divert(-1)dnl


#------------------------------------------------------------------------------

# Additional safety measures

undefine([substr])dnl   # Avoid collision with C/C++ naming

#------------------------------------------------------------------------------
divert(0)dnl
