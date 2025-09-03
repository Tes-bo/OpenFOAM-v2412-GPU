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
#     Logic rules, using bool or Foam::scalar for its storage.
#
#     `rule_logical_and'
#     `rule_logical_or'
#     `rule_logical_negate'
#     `rule_rule_ternary_op'
#     `rule_binary_logical_op'
#     `rule_mag_logical'
#
# Defined after inclusion
#      `_new_NNNfield'
#
# Defined here, but can be redefined later
#      `_logic_type_' `_logic_true_' `_logic_false_'
#
#------------------------------------------------------------------------------

# If previously defined, remove now
undefine([_logic_type_])dnl     # Eg, bool, Foam::scalar
undefine([_logic_true_])dnl     # Eg, true, Foam::scalar(1)
undefine([_logic_false_])dnl    # Eg, false, Foam::scalar(0)

#------------------------------------------------------------------------------
# rule_logical_and(inOut, [logicType])
# rule_logical_or(inOut, [logicType])
# rule_logical_negate(inOut, [logicType])
#
# Description
#     Production rule for basic logical operations
#
# Example - using volScalarField for selections
#
# rule_logical_and(lfield, Foam::scalar)
#
# lfield (lhs) ::= lfield (a) LAND lfield (b) .
# {
#     lhs = a;
#     Foam::FieldOps::assign
#     (
#         *lhs, (*a), make_obj(b),
#         Foam::expressions::boolAndOp<Foam::scalar>()
#     );
# }
#
# Example - using boolField for selections
#
# rule_logical_or(lfield) or rule_logical_or(lfield, bool)
#
# vfield (lhs) ::= lfield(a) LOR lfield (b) .
# {
#     lhs = a;
#
#     Foam::FieldOps::assign
#     (
#         *lhs, make_obj(cond), *a, make_obj(b),
#         Foam::orOp<bool>()
#     );
# }
#------------------------------------------------------------------------------

define([rule_logical_and],
[$1 (lhs) ::= $1 (a) LAND $1 (b) .
{
    lhs = a;
    Foam::FieldOps::assign
    (
        *lhs, *a, make_obj(b),
        Foam::expressions::boolAndOp<$2>()
    );
}]
)

define([rule_logical_or],
[$1 (lhs) ::= $1 (a) LOR $1 (b) .
{
    lhs = a;
    Foam::FieldOps::assign
    (
        *lhs, *a, make_obj(b),
        Foam::expressions::boolOrOp<$2>()
    );
}]
)

define([rule_logical_negate],
[$1 (lhs) ::= LNOT $1 (a). _lemon_precedence(NEGATE)
{
    lhs = a;
    Foam::FieldOps::assign
    (
        *lhs, *a,
        Foam::expressions::boolNotOp<$2>()
    );
}]
)


#------------------------------------------------------------------------------
# rule_ternary_op(inOut, cond, baseType)
#
# Description
#     Production rule for ternary field operations
#
# Uses [_logic_type_]
#
# Example - using boolField for selections
#
# rule_ternary_op(vfield, lfield, Foam::vector)
#
# vfield (lhs) ::= lfield(cond) QUESTION vfield (a) COLON vfield (b) .
# {
#     lhs = a;
#
#     Foam::FieldOps::ternarySelect<Foam::vector>
#     (
#         *lhs, make_obj(cond), *a, make_obj(b)
#     );
# }
#
# Example - using volScalarField for selections
#
# rule_ternary_op(vfield, lfield, Foam::vector)
#
# vfield (lhs) ::= lfield(cond) QUESTION vfield (a) COLON vfield (b) .
# {
#     lhs = a;
#
#     Foam::FieldOps::ternarySelect<Foam::vector, Foam::scalar>
#     (
#         *lhs, make_obj(cond), *a, make_obj(b),
#         Foam::expressions::boolOp<Foam::scalar>()
#     );
# }
#------------------------------------------------------------------------------

define([rule_ternary_op],
[$1 (lhs) ::= $2 (cond) QUESTION $1 (a) COLON $1 (b) .
{
    lhs = a;

    Foam::FieldOps::ternarySelect
    ifelse(_logic_type_, [bool],
    [<$3>(*lhs, make_obj(cond), *a, make_obj(b));],
    [<$3,_logic_type_>
    (
        *lhs, make_obj(cond), *a, make_obj(b),
        Foam::expressions::boolOp<_logic_type_>()
    );])
}]
)


#------------------------------------------------------------------------------
# rule_cast_logical(out, in, valType)
#
# Description
#     Production rules for type to logic functional casting
#
# Example
#     rule_cast_logical(lfield, sfield, Foam::scalar)
#------------------------------------------------------------------------------

define([rule_cast_logical],
[$1 (lhs) ::= BOOL LPAREN $2 (a) RPAREN .
{
ifelse($1, $2,
[    lhs = a;],
[    lhs = _new_$1();

    Foam::FieldOps::assign
    (
        *lhs, make_obj(a),
        Foam::expressions::boolOp<$3>()
    );])
}]
)


#------------------------------------------------------------------------------
# rule_mag_logical(out, in)
#
# Description
#     Production rules for logic to scalar
#
# Example
#     rule_mag_logical(sfield, lfield)
#------------------------------------------------------------------------------

define([rule_mag_logical],
[$1 (lhs) ::= MAG LPAREN $2 (a) RPAREN .
{
    lhs = _new_$1();
    Foam::FieldOps::assign
    (
        *lhs,
        make_obj(a),
        Foam::expressions::boolOp<_logic_type_>()
    );
}]
)


#------------------------------------------------------------------------------
# rule_binary_logical_op(out, in, valType, compare, tok, op)
#
# Description
#     Production rule for binary logical field operations
#     Operates on identical field types, producing a 'lfield'.
#
# Example
# rule_binary_logical_op(lfield, vfield, Foam::vector, greaterOp, GREATER, >)
#
# lfield (lhs) ::= vfield (a) GREATER vfield (b) .
# {
#    /// checkSizes("(a $6 b"), a, b);
#    lhs = driver->newVolField<Foam::scalar>().ptr();
#
#    Foam::FieldOps::assign
#    (*lhs, make_obj(a), make_obj(b), Foam::greaterOp<Foam::vector>());
# }
#------------------------------------------------------------------------------

define([rule_binary_logical_op],
[$1 (lhs) ::= $2 (a) $5 $2 (b) .
{
    lhs = _new_$1();
    Foam::FieldOps::assign(*lhs, make_obj(a), make_obj(b), Foam::$4<$3>());
}]
)


#------------------------------------------------------------------------------
# Logic constants
#
# _logic_type_, _logic_true_, _logic_false_
#------------------------------------------------------------------------------

define([use_bool_logic],
[define([_logic_type_], [bool])dnl
define([_logic_true_], [true])dnl
define([_logic_false_], [false])]
)

define([use_scalar_logic],
[define([_logic_type_], [Foam::scalar])dnl
define([_logic_true_], [Foam::scalar(1)])dnl
define([_logic_false_], [Foam::scalar(0)])]
)

# Define now as using Foam::scalar for the logic
use_scalar_logic()

#------------------------------------------------------------------------------
divert(0)dnl
