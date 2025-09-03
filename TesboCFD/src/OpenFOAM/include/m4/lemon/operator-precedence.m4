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
#     Defines standard operator precedence macro for lemon grammar.
#     Follows https://en.cppreference.com/w/cpp/language/operator_precedence
#
#------------------------------------------------------------------------------

define([operator_precedence],
[// [https://en.cppreference.com/w/cpp/language/operator_precedence]
%right QUESTION COLON .                 // 13: right-to-left
%left LOR  .                            // 12:
%left LAND .                            // 11:
%left BIT_OR  .                         // 10 (unused)
%left BIT_XOR .                         // 9
%left BIT_AND .                         // 8
%left EQUAL NOT_EQUAL .                 // 7
%left LESS LESS_EQUAL GREATER GREATER_EQUAL . // 6
%left PLUS MINUS .                      // 4
%left TIMES DIVIDE PERCENT .            // 3
%right NEGATE LNOT BIT_NOT .            // 2: right-to-left
%left DOT  .                            // 1: (method)]
)


define([standard_tokens],
[// Standard tokens for operators, constants and common types]
%token
  LPAREN RPAREN COMMA
  QUESTION COLON LOR LAND LNOT
  BIT_OR BIT_XOR BIT_AND BIT_NOT
  EQUAL NOT_EQUAL
  LESS LESS_EQUAL GREATER GREATER_EQUAL
  PLUS MINUS TIMES DIVIDE PERCENT
  NEGATE DOT
  BOOL LTRUE LFALSE
  NUMBER ZERO IDENTIFIER
.
)

#------------------------------------------------------------------------------
divert(0)dnl
