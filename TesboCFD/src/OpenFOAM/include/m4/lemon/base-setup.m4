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
#     A collection of 'base' setup of m4 macros for lemon, and setup
#     of quoting - uses square brackets (as per autoconf).
#
#     Specifies some reserved names that are to be defined prior to
#     invoking macros. These reserved names correspond to our internal
#     canonical names for basic data or field types used in OpenFOAM.
#
# Macros
#     begin_c_comment, end_c_comment
#     _lemon_precedence
#     compiler_pragmas
#     tmp_management
#
# Names for constants
#     _logic_type_, _logic_true_, _logic_false_
#
# Names for operands
#     _logic_, _scalar_, _vector_,
#     _sphTensor_, _symTensor_, _tensor_
#
# Values for the currently targeted rule
#     _target_, _value_type_, _scalar_arg_
#
# Note
#     The `undefine' occur immediately upon inclusion of this file.
#     Macro quoting is changed by this file!
#
# Examples
#
#     `lfield', `sfield', for operands
#     `vfield', `Foam::vector' for the `target' and `value_type'
#
#------------------------------------------------------------------------------

changequote(`[',`]')dnl Bracket quoting (as per autoconf)


#------------------------------------------------------------------------------
# begin_c_comment()
# end_c_comment()
#
# Description
#     Text for begin/end of C comments.
#     Can be useful if the m4 comment has been redefined to use C/C++ comments.
#------------------------------------------------------------------------------

define([begin_c_comment], [/*])
define([end_c_comment], [*/])


#------------------------------------------------------------------------------
# _lemon_precedence(name)
#
# Description
#     Add lemon precedence rule suffix. Eg, [NEGATE]
#------------------------------------------------------------------------------

define([_lemon_precedence], [ifelse([$#], [0], [], [[[$1]]])])


#------------------------------------------------------------------------------
# compiler_pragmas()
#
# Description
#     Common compiler pragmas when using lemon
#------------------------------------------------------------------------------

define([compiler_pragmas],
[dnl
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"])


#------------------------------------------------------------------------------
# tmp_management()
#
# Description
#     Defines file-local C++ functions to create a `Foam::tmp'
#     or a regular object from a pointer. This function wrapping
#     provides an effective means of automatically managing pointer
#     deletion.
#------------------------------------------------------------------------------

define([tmp_management],
[dnl
//- Create a tmp from a pointer, taking ownership
template<class T>
Foam::tmp<T> make_tmp(T* p)
{
    return Foam::tmp<T>(p);
}

//- Default [make_obj] is pass-through
template<class T>
const T& make_obj(const T& o) noexcept
{
    return o;
}

//- Move construct an object from a pointer and destroy the pointer
template<class T>
T make_obj(T*& p)
{
    T o(std::move(*p));
    delete p;
    p = nullptr;  // Prevent caller from deleting too
    return o;
}]
)


#------------------------------------------------------------------------------
# Clear reserved names
#
#------------------------------------------------------------------------------

# Logic

undefine([_logic_type_])
undefine([_logic_true_])
undefine([_logic_false_])

# Base types

undefine([_logic_])
undefine([_scalar_])
undefine([_vector_])
undefine([_sphTensor_])
undefine([_symTensor_])
undefine([_tensor_])

# Current target information

undefine([_target_])
undefine([_value_type_])
undefine([_scalar_arg_])

#------------------------------------------------------------------------------
divert(0)dnl
