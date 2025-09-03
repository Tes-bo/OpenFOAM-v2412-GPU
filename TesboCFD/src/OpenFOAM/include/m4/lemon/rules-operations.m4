divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019-2022 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Description
#     Collection of `standard' operations and type-specific ones.
#
#     `rules_compare_operations'
#     `rules_logical_operations'
#     `rules_scalar_operations'
#     `rules_vector_operations'
#     `rules_sphTensor_operations'
#     `rules_symTensor_operations'
#     `rules_tensor_operations'
#
# Defined after inclusion
#     `rule_binary_op'
#     `rule_ternary_op'
#     `rule_negate_op'
#     `rule_logical_and'
#     `rule_logical_or'
#     `rule_logical_negate'
#     `rule_binary_logical_op'
#
#------------------------------------------------------------------------------

# These are to be defined *after* inclusion - undefine now

undefine([rule_binary_op])
undefine([rule_ternary_op])
undefine([rule_logical_and])
undefine([rule_logical_or])
undefine([rule_logical_negate])
undefine([rule_binary_logical_op])
undefine([_logic_type_])


#------------------------------------------------------------------------------
# rules_compare_operations(out, in, valType)
#
# Description
#     Production rules for binary field operations, operating on
#     identical field types, producing a logical field.
#
# Example
# rules_compare_operations(lfield, vfield, Foam::vector)
#
# Depends on [_logic_type_]
#
#------------------------------------------------------------------------------

define([rules_compare_operations],
[rule_binary_logical_op($1, $2, $3, equalOp, EQUAL, ==)
rule_binary_logical_op($1, $2, $3, notEqualOp, NOT_EQUAL, !=)
rule_binary_logical_op($1, $2, $3, lessOp, LESS, <)
rule_binary_logical_op($1, $2, $3, lessEqualOp, LESS_EQUAL, <=)
rule_binary_logical_op($1, $2, $3, greaterOp, GREATER, >)
rule_binary_logical_op($1, $2, $3, greaterEqualOp, GREATER_EQUAL, >=)
rule_ternary_op($2, $1, $3)
])


#------------------------------------------------------------------------------
# rules_logical_operations(inOut, logicType)
#
# Description
#     Production rules for logical field operations.
#     Negation ('!') may need special treatment if Foam::scalar is used
#     for storage. Eg, treat zero-like (mag < 0.5) as true.
#
# Example
# rules_logical_operations(lfield, Foam::scalar)
#
#------------------------------------------------------------------------------

define([rules_logical_operations],
[$1 (lhs) ::= LPAREN $1 (a) RPAREN. { lhs = a; }
rule_logical_negate($@)
rule_logical_and($@)
rule_logical_or($@)
])


#------------------------------------------------------------------------------
# Operations - scalar
#
# Uses: _scalar_
#------------------------------------------------------------------------------

define([rules_scalar_operations],
[dnl
rule_binary_op(_scalar_, _scalar_, _scalar_, PLUS, +)
rule_binary_op(_scalar_, _scalar_, _scalar_, MINUS, -)
rule_binary_op(_scalar_, _scalar_, _scalar_, TIMES, *)
dnl
rule_scalar_divide(_scalar_, _scalar_, _scalar_)
rule_scalar_modulo(_scalar_)
dnl
dnl VectorSpace operations
dnl
rule_binary_op(_scalar_, _vector_, _vector_, BIT_AND, &)dnl> [inner]
])

#------------------------------------------------------------------------------
# Operations - vector
#
# Uses: _scalar_, _vector_
#------------------------------------------------------------------------------

define([rules_vector_operations],
[dnl
rule_binary_op(_vector_, _vector_, _vector_, PLUS, +)
rule_binary_op(_vector_, _vector_, _vector_, MINUS, -)
rule_binary_op(_vector_, _vector_, _scalar_, TIMES, *)
rule_binary_op(_vector_, _scalar_, _vector_, TIMES, *)
dnl
rule_scalar_divide(_vector_, _vector_, _scalar_)
dnl
dnl VectorSpace operations
dnl
rule_binary_op(_vector_, _vector_, _vector_, BIT_XOR, ^)
dnl
rule_binary_op(_vector_, _vector_, _tensor_, BIT_AND, &)
rule_binary_op(_vector_, _vector_, _symTensor_, BIT_AND, &)
rule_binary_op(_vector_, _vector_, _sphTensor_, BIT_AND, &)
rule_binary_op(_vector_, _tensor_, _vector_, BIT_AND, &)
rule_binary_op(_vector_, _symTensor_, _vector_, BIT_AND, &)
rule_binary_op(_vector_, _sphTensor_, _vector_, BIT_AND, &)
])

#------------------------------------------------------------------------------
# Operations - sphericalTensor
#
# Uses: _scalar_, _tensor_, _symTensor_, _sphTensor_
#------------------------------------------------------------------------------

define([rules_sphTensor_operations],
[dnl
rule_binary_op(_sphTensor_, _sphTensor_, _sphTensor_, PLUS, +)
rule_binary_op(_sphTensor_, _sphTensor_, _sphTensor_, MINUS, -)
rule_binary_op(_sphTensor_, _sphTensor_, _scalar_, TIMES, *)
rule_binary_op(_sphTensor_, _scalar_, _sphTensor_, TIMES, *)
dnl
rule_scalar_divide(_sphTensor_, _sphTensor_, _scalar_)
dnl
dnl VectorSpace operations
dnl
dnl TODO? (double inner):
dnl     rule_binary_op(_sphTensor_, _sphTensor_, _sphTensor_, LAND, &&)
])

#------------------------------------------------------------------------------
# Operations - symmTensor
#
# Uses: _scalar_, _tensor_, _symTensor_, _sphTensor_
#------------------------------------------------------------------------------

define([rules_symTensor_operations],
[dnl
rule_binary_op(_symTensor_, _symTensor_, _symTensor_, PLUS, +)
rule_binary_op(_symTensor_, _symTensor_, _sphTensor_, PLUS, +)
rule_binary_op(_symTensor_, _sphTensor_, _symTensor_, PLUS, +)
dnl
rule_binary_op(_symTensor_, _symTensor_, _symTensor_, MINUS, -)
rule_binary_op(_symTensor_, _symTensor_, _sphTensor_, MINUS, -)
rule_binary_op(_symTensor_, _sphTensor_, _symTensor_, MINUS, -)
dnl
rule_binary_op(_symTensor_, _symTensor_, _scalar_, TIMES, *)
rule_binary_op(_symTensor_, _scalar_, _symTensor_, TIMES, *)
dnl
rule_scalar_divide(_symTensor_, _symTensor_, _scalar_)
dnl
dnl VectorSpace operations
dnl
rule_binary_op(_symTensor_, _symTensor_, _sphTensor_, BIT_AND, &)
rule_binary_op(_symTensor_, _sphTensor_, _symTensor_, BIT_AND, &)
])

#------------------------------------------------------------------------------
# Operations - tensor
#
# Uses: _scalar_, _tensor_, _symTensor_, _sphTensor_
#------------------------------------------------------------------------------

define([rules_tensor_operations],
[dnl
rule_binary_op(_tensor_, _tensor_, _tensor_, PLUS, +)
rule_binary_op(_tensor_, _tensor_, _symTensor_, PLUS, +)
rule_binary_op(_tensor_, _symTensor_, _tensor_, PLUS, +)
rule_binary_op(_tensor_, _tensor_, _sphTensor_, PLUS, +)
rule_binary_op(_tensor_, _sphTensor_, _tensor_, PLUS, +)
dnl
rule_binary_op(_tensor_, _tensor_, _tensor_, MINUS, -)
rule_binary_op(_tensor_, _tensor_, _symTensor_, MINUS, -)
rule_binary_op(_tensor_, _symTensor_, _tensor_, MINUS, -)
rule_binary_op(_tensor_, _tensor_, _sphTensor_, MINUS, -)
rule_binary_op(_tensor_, _sphTensor_, _tensor_, MINUS, -)
dnl
rule_binary_op(_tensor_, _tensor_, _scalar_, TIMES, *)
rule_binary_op(_tensor_, _scalar_, _tensor_, TIMES, *)
dnl
rule_scalar_divide(_tensor_, _tensor_, _scalar_)
dnl
dnl VectorSpace operations
dnl
rule_binary_op(_tensor_, _vector_, _vector_, TIMES, *)
dnl
rule_binary_op(_tensor_, _tensor_, _tensor_, BIT_AND, &)
rule_binary_op(_tensor_, _tensor_, _sphTensor_, BIT_AND, &)
rule_binary_op(_tensor_, _tensor_, _symTensor_, BIT_AND, &)
rule_binary_op(_tensor_, _sphTensor_, _tensor_, BIT_AND, &)
rule_binary_op(_tensor_, _symTensor_, _tensor_, BIT_AND, &)
rule_binary_op(_tensor_, _symTensor_, _symTensor_, BIT_AND, &)
dnl
])


#------------------------------------------------------------------------------
divert(0)dnl
