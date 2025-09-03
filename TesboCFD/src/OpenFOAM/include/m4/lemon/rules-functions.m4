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
#     Collection of `standard' functions and type-specific ones.
#
#     `rules_inplace_unary'
#     `rules_scalar_functions'
#     `rules_vector_functions'
#     `rules_sphTensor_functions'
#     `rules_symTensor_functions'
#     `rules_tensor_functions'
#
# Defined after inclusion
#     `rule_unary_func'         `incomplete_rule_unary_func'
#     `rule_binary_func'
#     `rule_inplace_unary'
#
#------------------------------------------------------------------------------

# These are to be defined *after* inclusion - undefine now

undefine([rule_unary_func])
undefine([rule_binary_func])
undefine([rule_inplace_unary])
undefine([incomplete_rule_unary_func])


#------------------------------------------------------------------------------
# rules_inplace_unary(inOut)
#
# Description
#     Production rules for standard unary field reductions.
#
#     Calls [rule_inplace_unary()] with
#     (MIN, min), (MAX, max), (SUM, sum), (AVERAGE, average)
#
# Example
#     rules_inplace_unary(sfield)
#------------------------------------------------------------------------------

define([rules_inplace_unary],
[rule_inplace_unary($1, MIN, Foam::min)
rule_inplace_unary($1, MAX, Foam::max)
rule_inplace_unary($1, SUM, Foam::sum)
rule_inplace_unary($1, AVERAGE, Foam::average)]
)


#------------------------------------------------------------------------------
# rules_inplace_gUnary(inOut)
#
# Description
#     Production rules for standard global unary field reductions.
#
#     Calls [rule_inplace_unary()] with
#     (MIN, gMin), (MAX, gMax), (SUM, gSum), (AVERAGE, gAverage)
#
# Example
#     rules_inplace_gUnary(sfield)
#------------------------------------------------------------------------------

define([rules_inplace_gUnary],
[rule_inplace_unary($1, MIN, Foam::gMin)
rule_inplace_unary($1, MAX, Foam::gMax)
rule_inplace_unary($1, SUM, Foam::gSum)
rule_inplace_unary($1, AVERAGE, Foam::gAverage)]
)


#------------------------------------------------------------------------------
# Functions - magnitude
#------------------------------------------------------------------------------

define([rules_mag_functions],
[dnl
rule_unary_func($1, $2, MAG, Foam::mag)
rule_unary_func($1, $2, MAGSQR, Foam::magSqr)
])


#------------------------------------------------------------------------------
# Functions - scalar
#------------------------------------------------------------------------------

define([rules_scalar_functions],
[dnl
rule_unary_func(_scalar_, _scalar_, EXP, Foam::exp)
rule_unary_func(_scalar_, _scalar_, LOG, Foam::log)
rule_unary_func(_scalar_, _scalar_, LOG10, Foam::log10)
rule_unary_func(_scalar_, _scalar_, SQR, Foam::sqr)
rule_unary_func(_scalar_, _scalar_, SQRT, Foam::sqrt)
rule_unary_func(_scalar_, _scalar_, CBRT, Foam::cbrt)
rule_unary_func(_scalar_, _scalar_, SIN, Foam::sin)
rule_unary_func(_scalar_, _scalar_, COS, Foam::cos)
rule_unary_func(_scalar_, _scalar_, TAN, Foam::tan)
rule_unary_func(_scalar_, _scalar_, ASIN, Foam::asin)
rule_unary_func(_scalar_, _scalar_, ACOS, Foam::acos)
rule_unary_func(_scalar_, _scalar_, ATAN, Foam::atan)
rule_unary_func(_scalar_, _scalar_, SINH, Foam::sinh)
rule_unary_func(_scalar_, _scalar_, COSH, Foam::cosh)
rule_unary_func(_scalar_, _scalar_, TANH, Foam::tanh)
rule_binary_func(_scalar_, _scalar_, _scalar_, POW, Foam::pow)
rule_binary_func(_scalar_, _scalar_, _scalar_, ATAN2, Foam::atan2)
dnl
incomplete_rule_unary_func(_scalar_, _scalar_, POS, Foam::pos)
incomplete_rule_unary_func(_scalar_, _scalar_, NEG, Foam::neg)
incomplete_rule_unary_func(_scalar_, _scalar_, POS0, Foam::pos0)
incomplete_rule_unary_func(_scalar_, _scalar_, NEG0, Foam::neg0)
incomplete_rule_unary_func(_scalar_, _scalar_, SIGN, Foam::sign)
])

#------------------------------------------------------------------------------
# Functions - vector
#------------------------------------------------------------------------------

define([rules_vector_functions], [])


#------------------------------------------------------------------------------
# Functions - sphericalTensor
#------------------------------------------------------------------------------

define([rules_sphTensor_functions], [])


#------------------------------------------------------------------------------
# Functions - symmTensor
#------------------------------------------------------------------------------

define([rules_symTensor_functions],[])


#------------------------------------------------------------------------------
# Functions - tensor
#------------------------------------------------------------------------------

define([rules_tensor_functions], [])


#------------------------------------------------------------------------------
divert(0)dnl
