divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Description
#     Collection of common functions that should work on all (non-logical)
#     field types
#
#     `rules_standard'
#
# Defined after inclusion
#     `rules_compare_operations'
#     `rule_negate_op'
#     `rule_const_multiply'
#     `rule_binary_func'
#
#------------------------------------------------------------------------------

# These are to be defined *after* inclusion - undefine now

undefine([rules_compare_operations])
undefine([rule_negate_op])
undefine([rule_const_multiply])
undefine([rule_binary_func])


#------------------------------------------------------------------------------
# rules_standard(target, valType, logicFieldType)
#
# Description
#     Production rules for some common functions that should work on
#     all (non-logical) field types
#
# Uses
# - [rule_negate_op]
# - [rule_binary_func()] with  (MIN, min), (MAX, max)
# - [rule_const_multiply()] for degToRad, radToDeg
# - [rules_compare_operations]
#
# Example
#     rules_standard(sfield, Foam::scalar, lfield)
#------------------------------------------------------------------------------

define([rules_standard],
[dnl
$1 (lhs) ::= LPAREN $1 (a) RPAREN. { lhs = a; }
rule_negate_op($@)
dnl
rule_const_multiply($1, $1, Foam::degToRad(), DEG_TO_RAD)
rule_const_multiply($1, $1, Foam::radToDeg(), RAD_TO_DEG)
rule_binary_func($1, $1, $1, MIN, Foam::min)
rule_binary_func($1, $1, $1, MAX, Foam::max)
dnl
rules_compare_operations($3, $1, $2)
]
)


#---------------------------------------------------------------------------
divert(0)dnl
