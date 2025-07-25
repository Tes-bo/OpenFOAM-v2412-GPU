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
#     Various "boilerplate" parser methods (C++)
#
# Requires
#     IOmanip.H
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# parser_code_static_methods(clsName)
#
# Description
#     Emit common parser static methods
#        `tokenName`
#        `printTokenNames`
#        `printRules`
#
# Note
#     Array access uses `*(array + i)` to avoid [] square brackets,
#     which may be used for m4 quoting, unless we switched back to `' !!
#
# Example
# parser_code_static_methods(Foam::expressions::fieldExpr::parser)
#
#------------------------------------------------------------------------------

define([parser_code_static_methods],
[dnl
Foam::word $1::tokenName(int i)
{
    #ifndef NDEBUG
    if (i > 0 && unsigned(i) < (sizeof(yyTokenName) / sizeof(char*)))
    {
        return *(yyTokenName + i);
    }
    return "<invalid>";
    #else
    return "";
    #endif
}

void $1::printTokenNames(Ostream& os)
{
    #ifndef NDEBUG
    const unsigned nElem(sizeof(yyTokenName) / sizeof(char*));
    for (unsigned i = 1; i < nElem; ++i) // start = 1 (skip termination token)
    {
        os << *(yyTokenName + i) << nl;
    }
    #endif
}

void $1::printRules(Ostream& os)
{
    #ifndef NDEBUG
    const unsigned nElem(sizeof(yyRuleName) / sizeof(char*));

    // Easy way to count number of digits
    const unsigned width(std::to_string(nElem).length());

    for (unsigned i = 0; i < nElem; ++i)
    {
        os << setw(width) << i << ": " << *(yyRuleName + i) << nl;
    }
    #endif
}]
)

#------------------------------------------------------------------------------
divert(0)dnl
