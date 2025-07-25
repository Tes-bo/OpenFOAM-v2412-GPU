divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019-2021 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Description
#     Field handling m4/lemon macros.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# declare_field(name, fieldType, value_type, newMethod, getMethod)
#
# Description
#     Define lemon parser field types and macros
#         _new_NNfield
#         _get_NNfield
#         _value_type_NNfield
#
#     The new and get macros shall deliver the correct type
#     (eg, pointer) for use with lemon.
#
# Example
#     declare_field
#     (
#         sfield,
#         Foam::volScalarField,
#         Foam::scalar,
#         newVolField,
#         getVolField
#     )
#
# %type sfield { Foam::volScalarField* }
#
# def: _new_sfield
# def: _get_sfield
# def: _value_type_sfield
#------------------------------------------------------------------------------

define([declare_field],
[dnl
define([_value_type_]$1, [$3])dnl
define([_new_]$1, [driver->$4<$3>($][1).ptr()])dnl
define([_get_]$1, [driver->$5<$3>($][1).ptr()])dnl
%type $1 { $2* }])


#------------------------------------------------------------------------------
# rule_get_field(out, tok)
#
# Description
#     Production rule for driver get method
#
# Example
# rule_get_field(sfield, SCALAR_ID)
#
# sfield (lhs) ::= SCALAR_ID (name) .
# {
#     lhs = driver->getVolField<Foam::scalar>(make_obj(name.name_)).ptr();
# }
#------------------------------------------------------------------------------

define([rule_get_field],
[$1 (lhs) ::= $2 (name) .
{
    lhs = _get_$1(make_obj(name.name_));
}]
)


#------------------------------------------------------------------------------
# rule_field_from_value(out, in, [prefix])
#
# Description
#     Production rule for value to field promotion
#
# Example
# rule_field_from_value(sfield, svalue)
# rule_field_from_value(psfield, svalue, POINT_EXPR)
#
# sfield (lhs) ::= svalue (a) .
# {
#     lhs = driver->newField<Foam::scalar>(make_obj(a)).ptr();
# }
#
# psfield (lhs) ::= POINT_EXPR LPAREN svalue (a) RPAREN .
# {
#     lhs = driver->newPointField<Foam::scalar>(make_obj(a)).ptr();
# }
#------------------------------------------------------------------------------

define([rule_field_from_value],
[$1 (lhs) ::= ifelse($3,[],[],[$3 LPAREN ])$2 (a)ifelse($3,[],[],[ RPAREN]) .
{
    lhs = _new_$1(make_obj(a));
}]
)


#------------------------------------------------------------------------------
# rule_negate_op(target)
#
# Description
#     Production rules for field negation
#
# Example
#     rule_negate_op(sfield)
#------------------------------------------------------------------------------

define([rule_negate_op],
[$1 (lhs) ::= MINUS $1 (a) . _lemon_precedence(NEGATE)
{
    lhs = a; lhs->negate();
}]
)


#------------------------------------------------------------------------------
# rule_binary_op(out, in1, in2, tok, op)
#
# Description
#     Production rule for binary field operations
#
# Example
# rule_binary_op(sfield, sfield, sfield, PLUS, +)
#
# sfield (lhs) ::= sfield (a) PLUS sfield (b) .
# {
#     /// checkSizes("(a $5 b)", a, b);
#     lhs = (make_tmp(a) + make_tmp(b)).ptr();
# }
#
#------------------------------------------------------------------------------

define([rule_binary_op],
[$1 (lhs) ::= $2 (a) $4 $3 (b) .
{
    lhs = (make_tmp(a) $5 make_tmp(b)).ptr();
}]
)


#------------------------------------------------------------------------------
# rule_inplace_unary(inOut, tok, func)
#
# Description
#     Production rule for inplace unary field reductions
#
# Example
# rule_inplace_unary(vfield, MIN, Foam::min)
#
# vfield (lhs) ::= MIN LPAREN vfield (a) RPAREN .
# {
#     lhs = a; *lhs = Foam::min (*lhs);
# }
#------------------------------------------------------------------------------

define([rule_inplace_unary],
[$1 (lhs) ::= $2 LPAREN $1 (a) RPAREN .
{
    lhs = a; *lhs = $3 (*lhs);
}]
)

#------------------------------------------------------------------------------
# rule_unary_func(out, in1, tok, func)
#
# Description
#     Production rule for unary functions
#
# Example
# rule_unary_func(sfield, vfield, MAG, Foam::mag)
#
# sfield (lhs) ::= MAG LPAREN vfield (a) RPAREN .
# {
#     lhs = Foam::mag(make_tmp(a)).ptr();
# }
#------------------------------------------------------------------------------

define([rule_unary_func],
[$1 (lhs) ::= $3 LPAREN $2 (a) RPAREN .
{
    lhs = $4 (make_tmp(a)).ptr();
}]
)


#------------------------------------------------------------------------------
# rule_binary_func(out, in1, in2, tok, func)
#
# Description
#     Production rule for binary functions
#
# Example
# rule_binary_func(sfield, sfield, sfield, POW, Foam::pow)
#
# sfield (lhs) ::= POW LPAREN sfield (a) COMMA sfield (b) RPAREN .
# {
#     /// checkSizes("$5 (a, b)", a, b);
#     lhs = Foam::pow(make_tmp(a), make_tmp(b)).ptr();
# }
#------------------------------------------------------------------------------

define([rule_binary_func],
[$1 (lhs) ::= $4 LPAREN $2 (a) COMMA $3 (b) RPAREN .
{
    lhs = $5(make_tmp(a), make_tmp(b)).ptr();
}]
)


#------------------------------------------------------------------------------
# rule_const_multiply(out, in, multiplier, tok)
#
# Description
#     Production rule for constant multiplier method
#
# Example
# rule_const_multiply(sfield, sfield, Foam::degToRad(), DEG_TO_RAD)
#
# sfield(lhs) ::= DEG_TO_RAD LPAREN sfield(a) RPAREN .
# {
#     lhs = ((Foam::degToRad()) * make_tmp(a)).ptr();
# }
#------------------------------------------------------------------------------

define([rule_const_multiply],
[$1 (lhs) ::= $4 LPAREN $2 (a) RPAREN .
{
    lhs = (($3) * make_tmp(a)).ptr();
}]
)


#------------------------------------------------------------------------------
# rule_scalar_divide(out, in1, in2, [value_type])
#
# Description
#     Production rule for division by scalar operation
#
# Example
# rule_scalar_divide(vfield, vfield, sfield, Foam::vector)
#
# vfield(lhs) ::= vfield (a) DIVIDE sfield (b) .
# {
#     lhs = _new_$1();
#     Foam::FieldOps::assign
#     (
#         *lhs,
#         make_obj(a),
#         make_obj(b),
#         Foam::scalarDivideOp<Foam::vector>()
#     );
# }
#------------------------------------------------------------------------------

define([rule_scalar_divide],
[$1 (lhs) ::= $2 (a) DIVIDE $3 (b) .
{
    lhs = _new_$1();
    Foam::FieldOps::assign
    (
        *lhs,
        make_obj(a),
        make_obj(b),
        Foam::scalarDivideOp<ifelse($4,[], _value_type_$1, $4)>()
    );
}]
)


#------------------------------------------------------------------------------
# rule_scalar_modulo(inOut)
#
# Description
#     Production rule for scalar modulo operator
#
# Example
#     rule_scalar_modulo(sfield)
#
# NOTE
#     Swak uses the equivalent of
#     Foam::expressions::exprDriverOps::swakModuloOp()
#------------------------------------------------------------------------------

define([rule_scalar_modulo],
[$1 (lhs) ::= $1 (a) PERCENT $1 (b) .
{
    lhs = _new_$1();
    Foam::FieldOps::assign
    (
        *lhs,
        make_obj(a),
        make_obj(b),
        Foam::scalarModuloOp<Foam::scalar>()
    );
}]
)


#------------------------------------------------------------------------------
# rule_unary_assign(out, in, tok, function)
#
# Description
#     Production rule for a unary function,
#     using FieldOps::assign for the implementation.
#
# Example
# rule_unary_assign(sfield, sfield, FLOOR, Foam::floorOp<Foam::scalar>())
#
#------------------------------------------------------------------------------

define([rule_unary_assign],
[$1 (lhs) ::= $3 LPAREN $2 (a) RPAREN .
{
    lhs = _new_$1();
    Foam::FieldOps::assign(*lhs, make_obj(a), $4);
}]
)


#------------------------------------------------------------------------------
# rule_binary_assign(out, in1, in2, tok, function)
#
# Description
#     Production rule for a binary function,
#     using FieldOps::assign for the implementation.
#
# Example
# rule_binary_assign(sfield, sfield, sfield, HYPOT, Foam::hypot)
#
#------------------------------------------------------------------------------

define([rule_binary_assign],
[$1 (lhs) ::= $4 LPAREN $2 (a) COMMA $3 (b) RPAREN .
{
    lhs = _new_$1();
    Foam::FieldOps::assign(*lhs, make_obj(a), make_obj(b), $5);
}]
)


#------------------------------------------------------------------------------
# rule_driver_nullary(out, tok, func)
#
# Description
#     Production rule for driver-specific nullary method.
#
# Example
# rule_driver_nullary(vfield, POS, field_cellCentre)
#
# vfield (lhs) ::= POS LPAREN RPAREN .
# {
#     lhs = driver->field_cellCentre().ptr();
# }
#------------------------------------------------------------------------------
define([rule_driver_nullary],
[$1 (lhs) ::= $2 LPAREN RPAREN .
{
    lhs = driver->$3().ptr();
}]
)


#------------------------------------------------------------------------------
# rule_driver_inplace_unary(inOut, tok, method)
#
# Description
#     Production rule for a driver-specific unary method
#     modifying the field inplace.
#
# Example
# rule_driver_inplace_unary(sfield, WEIGHTED_AVERAGE, volAverage)
#
# sfield(lhs) ::= WEIGHTED_AVERAGE LPAREN sfield (a) RPAREN .
# {
#     lhs = a; *lhs = driver->volAverage(*lhs);
# }
#------------------------------------------------------------------------------

define([rule_driver_inplace_unary],
[$1 (lhs) ::= $2 LPAREN $1 (a) RPAREN .
{
    lhs = a; *lhs = driver->$3(*lhs);
}]
)


#------------------------------------------------------------------------------
# rule_driver_unary(out, in, tok, method, [value_type])
#
# Description
#     Production rule for a driver-specific unary method
#
# Example
# rule_driver_unary(sfield, psfield, POINT_TO_CELL, pointToCell)
#
# sfield(lhs) ::= POINT_TO_CELL LPAREN psfield (a) RPAREN .
# {
#     lhs = driver->pointToCell(make_obj(a)).ptr();
# }
#------------------------------------------------------------------------------

define([rule_driver_unary],
[$1 (lhs) ::= $3 LPAREN $2 (a) RPAREN .
{
    lhs = driver->$4[]dnl       # The method call
ifelse($5,[],[],[<$5>])dnl      # Optional template parameter (value_type)
(make_obj(a)).ptr();
}]
)


#------------------------------------------------------------------------------
# rule_driver_unary_named(out, tok, identType, method, [value_type])
#
# Description
#     Production rule for a driver-specific unary method
#
# Example
#     rule_driver_unary_named
#     (
#        sfield,
#        SN_GRAD,
#        SCALAR_ID,
#        patchNormalField,
#        Foam::scalar
#    )
#
# sfield(lhs) ::= SN_GRAD LPAREN SCALAR_ID (name) RPAREN .
# {
#     lhs = driver->patchNormalField<Foam::scalar>(make_obj(name.name_)).ptr();
# }
#
#------------------------------------------------------------------------------

define([rule_driver_unary_named],
[$1 (lhs) ::= $2 LPAREN $3 (name) RPAREN .
{
    lhs = driver->$4[]dnl       # The method call
ifelse($5,[],[],[<$5>])dnl      # Optional template parameter (value_type)
(make_obj(name.name_)).ptr();
}]
)


#------------------------------------------------------------------------------
# Forward "incomplete" versions to allow overload (eg, point fields)

define([incomplete_rule_unary_func], defn([rule_unary_func]))
define([incomplete_rule_binary_func], defn([rule_binary_func]))


#------------------------------------------------------------------------------
divert(0)dnl
