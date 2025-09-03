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
#     Rules for vector/tensor `zip' functions, which are used to combine
#     scalar fields into vectors or tensors.
#
#     `rule_method_component'
#     `rule_tensor_transpose'
#     `rule_sphTensor_transpose'
#     `rule_symTensor_transpose'
#
#     `rule_vector_zip'
#     `rule_tensor_zip'
#     `rule_symTensor_zip'
#     `rule_sphTensor_zip'
#
#     `rule_tensor_unzipCol'
#     `rule_tensor_unzipRow'
#     `rule_tensor_unzipDiag'
#
#     Default `field_read_access' `field_write_access'
#
#     For pointFields the access macros can be redefined to use
#     the primitive fields.
#
# After inclusion - must define corresponding `_new_NNNfield' macros
#
#------------------------------------------------------------------------------

# These are to be defined *after* inclusion - undefine now

undefine([field_read_access])
undefine([field_write_access])

#------------------------------------------------------------------------------
# rule_method_component(out, in, tok, cmpt)
#
# Description
#     Production rule for method component
#
# Example
# rule_method_component(sfield, tfield, CMPT_XX, Foam::tensor::XX)
#
# sfield (lhs) ::= tfield (a) DOT method_call CMPT_XX LPAREN RPAREN. [DOT]
# {
#     lhs = make_tmp(a)->component(Foam::tensor::XX).ptr();
# }
#------------------------------------------------------------------------------

define([rule_method_component],
[$1 (lhs) ::= $2 (a) DOT $3 LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = make_tmp(a)->component($4).ptr();
}]
)


#------------------------------------------------------------------------------
# rule_sphTensor_transpose(source)
#
# Description
#     Transpose field
#------------------------------------------------------------------------------

define([rule_sphTensor_transpose],
[$1 (lhs) ::= $1 (a) DOT TRANSPOSE LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = a; // no-op
}]
)


#------------------------------------------------------------------------------
# rule_symTensor_transpose(source)
#
# Description
#     Transpose field
#------------------------------------------------------------------------------

define([rule_symTensor_transpose],
[$1 (lhs) ::= $1 (a) DOT TRANSPOSE LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = a; // no-op
}]
)


#------------------------------------------------------------------------------
# rule_tensor_transpose(source)
#
# Description
#     Transpose field
#------------------------------------------------------------------------------

define([rule_tensor_transpose],
[$1 (lhs) ::= $1 (a) DOT TRANSPOSE LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = a;
    Foam::T(*lhs, *a);
}]
)


#------------------------------------------------------------------------------
# rule_vector_zip(target, sources, tok)
#
# Description
#     Combine three scalars to produce a vector
#     - vector(x, y, z)
#
# A rule such as [x y z] does not reduce well - needs more investigation.
#
# Example
# rule_vector_zip(vfield, sfield, VECTOR)
#------------------------------------------------------------------------------

define([rule_vector_zip],
[$1 (lhs) ::= $3 LPAREN $2 (x) COMMA $2 (y) COMMA $2 (z) RPAREN.
{
    lhs = _new_$1();

    Foam::zip
    (
        field_write_access(*lhs),
        field_read_access(make_obj(x)),
        field_read_access(make_obj(y)),
        field_read_access(make_obj(z))
    );
}]
)


#------------------------------------------------------------------------------
# rule_sphTensor_zip(target, sources, tok)
#
# Description
#     Combine one scalar to produce a sphericalTensor
#     - sphericalTensor(ii)
#------------------------------------------------------------------------------

define([rule_sphTensor_zip],
[$1 (lhs) ::= $3 LPAREN $2 (ii) RPAREN.
{
    lhs = _new_$1();

    Foam::zip
    (
        field_write_access(*lhs),
        field_read_access(make_obj(ii))
    );
}]
)


#------------------------------------------------------------------------------
# rule_symTensor_zip(target, sources, tok)
#
# Description
#     Combine six scalars to produce a symmTensor
#     - symmTensor(xx, xy, xy, yy, yz, zz)
#------------------------------------------------------------------------------

define([rule_symTensor_zip],
[$1 (lhs) ::= $3 LPAREN
    $2 (xx) COMMA $2 (xy) COMMA $2 (xz) COMMA
    $2 (yy) COMMA $2 (yz) COMMA
    $2 (zz)
RPAREN.
{
    lhs = _new_$1();

    Foam::zip
    (
        field_write_access(*lhs),
        field_read_access(make_obj(xx)),
        field_read_access(make_obj(xy)),
        field_read_access(make_obj(xz)),
        field_read_access(make_obj(yy)),
        field_read_access(make_obj(yz)),
        field_read_access(make_obj(zz))
    );
}]
)


#------------------------------------------------------------------------------
# rule_tensor_zip(target, sources, tok)
#
# Description
#     Combine scalars to produce a tensor
#
# Example
# rule_tensor_zip(stfield, ssfield, TENSOR)
#------------------------------------------------------------------------------

define([rule_tensor_zip],
[$1 (lhs) ::= $3 LPAREN
    $2 (xx) COMMA $2 (xy) COMMA $2 (xz) COMMA
    $2 (yx) COMMA $2 (yy) COMMA $2 (yz) COMMA
    $2 (zx) COMMA $2 (zy) COMMA $2 (zz)
RPAREN.
{
    lhs = _new_$1();

    Foam::zip
    (
        field_write_access(*lhs),
        field_read_access(make_obj(xx)),
        field_read_access(make_obj(xy)),
        field_read_access(make_obj(xz)),
        field_read_access(make_obj(yx)),
        field_read_access(make_obj(yy)),
        field_read_access(make_obj(yz)),
        field_read_access(make_obj(zx)),
        field_read_access(make_obj(zy)),
        field_read_access(make_obj(zz))
    );
}]
)

#------------------------------------------------------------------------------
# rule_tensor_unzipDiag(out, in)
#
# Description
#     Extract vector diagonal from tensor, as a method call.
#
# Example
# rule_tensor_unzipDiag(vfield, tfield)
#
# vfield (lhs) ::= tfield (a) DOT DIAG LPAREN RPAREN. [DOT]
# {
#     lhs = _new_vfield();
#     Foam::unzipDiag(make_obj(a), *lhs);
# }
#------------------------------------------------------------------------------

define([rule_tensor_unzipDiag],
[$1 (lhs) ::= $2 (a) DOT DIAG LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = _new_$1();

    Foam::unzipDiag(field_read_access(make_obj(a)), field_write_access(*lhs));
}]
)


#------------------------------------------------------------------------------
# rule_tensor_unzipCol(out, in, tok, cmpt)
#
# Description
#     Extract vector column from tensor, as a method call.
#
# Example
# rule_tensor_unzipRow(vfield, tfield, CMPT_CX, Foam::vector::X)
#
# vfield (lhs) ::= tfield (a) DOT CMPT_CX LPAREN RPAREN. [DOT]
# {
#     lhs = _new_vfield();
#     Foam::unzipCol(make_obj(a), Foam::vector::X, *lhs);
# }
#------------------------------------------------------------------------------

define([rule_tensor_unzipCol],
[$1 (lhs) ::= $2 (a) DOT $3 LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = _new_$1();

    Foam::unzipCol
    (
        field_read_access(make_obj(a)),
        $4,
        field_write_access(*lhs)
    );
}]
)

#------------------------------------------------------------------------------
# rule_tensor_unzipRow(out, in, tok, cmpt)
#
# Description
#     Extract vector row from tensor, as a method call.
#
# Example
# rule_tensor_unzipRow(vfield, tfield, CMPT_X, Foam::vector::X)
#
# vfield (lhs) ::= tfield (a) DOT CMPT_X LPAREN RPAREN. [DOT]
# {
#     lhs = _new_vfield();
#     Foam::unzipRow(make_obj(a), Foam::vector::X, *lhs);
# }
#------------------------------------------------------------------------------

define([rule_tensor_unzipRow],
[$1 (lhs) ::= $2 (a) DOT $3 LPAREN RPAREN. _lemon_precedence(DOT)
{
    lhs = _new_$1();

    Foam::unzipRow
    (
        field_read_access(make_obj(a)),
        $4,
        field_write_access(*lhs)
    );
}]
)


#------------------------------------------------------------------------------
# Default field access rules
#------------------------------------------------------------------------------

define([field_write_access], [($1)])
define([field_read_access], [($1)])


#------------------------------------------------------------------------------
divert(0)dnl
