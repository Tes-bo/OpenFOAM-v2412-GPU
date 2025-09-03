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
#     Named versions (as m4 macros) of single quoted characters to avoid
#     quoting complexity when writing m4 macros of some bison rules.
#
#------------------------------------------------------------------------------

changequote(`[',`]')dnl Bracket quoting (as per autoconf)

#------------------------------------------------------------------------------

undefine([PLUS])define([PLUS], ['+'])
undefine([MINUS])define([MINUS], ['-'])
undefine([TIMES])define([TIMES], ['*'])
undefine([DIVIDE])define([DIVIDE], ['/'])
undefine([PERCENT])define([PERCENT], ['%'])dnl - eg modulo
undefine([QUESTION])define([QUESTION], ['?'])dnl - eg ternary
undefine([COLON])define([COLON], [':'])dnl - eg ternary
undefine([COMMA])define([COMMA], [','])

undefine([LPAREN])define([LPAREN], ['('])
undefine([RPAREN])define([RPAREN], [')'])
undefine([LBRACE])define([LBRACE], ['{'])
undefine([RBRACE])define([RBRACE], ['}'])


#------------------------------------------------------------------------------

changequote([`],['])dnl Standard m4 quoting

#------------------------------------------------------------------------------
divert(0)dnl
