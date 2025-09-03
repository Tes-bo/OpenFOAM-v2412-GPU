/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Flex-based parsing of STL ASCII format

\*---------------------------------------------------------------------------*/

%option prefix="yySTL"
%option yyclass="yySTLFlexLexer"

%{

 /* ------------------------------------------------------------------------ *\
   ------ local definitions
 \* ------------------------------------------------------------------------ */

#include "STLAsciiParse.H"
#include "STLReader.H"
#include "OSspecific.H"

// Flex may use register, which is deprecated and incompatible with C++17
#pragma clang diagnostic ignored "-Wdeprecated-register"

using namespace Foam;

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! \cond dummy
#if YY_FLEX_MAJOR_VERSION <= 2 && YY_FLEX_MINOR_VERSION <= 5 && YY_FLEX_SUBMINOR_VERSION < 34
int yyFlexLexer::yylex()
{
    FatalErrorInFunction
        << "Should not have called this function"
        << abort(FatalError);
    return 0;
}
#endif
//! \endcond

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! \cond dummy
#if YY_FLEX_MAJOR_VERSION <= 2 && YY_FLEX_MINOR_VERSION <= 5 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yySTLFlexLexer::yywrap()
#endif
{
    return 1;
}
//! \endcond


//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (zoneIds).
//  The facets are within a solid/endsolid grouping
class STLAsciiParseFlex
:
    public Detail::STLAsciiParse,
    public yySTLFlexLexer
{
    word startError_;

public:

    // Constructors

    //- From input stream, with the estimated number of triangles in the STL
    STLAsciiParseFlex(istream* is, const label nTrisEstimated)
    :
        Detail::STLAsciiParse(nTrisEstimated),
        yySTLFlexLexer(is)
    {}


    // Member Functions

    //- The lexer function itself
    int lex();

    //- Execute lexer
    void execute()
    {
        while (lex()) {}
    }
};


 /* ------------------------------------------------------------------------ *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------ */

#undef YY_DECL
#define YY_DECL int STLAsciiParseFlex::lex()

%}

white                 [ \t\f\r]
space                 {white}*
some_space            {white}+

alpha                 [_A-Za-z]
digit                 [0-9]

intNum                [-+]?{digit}+

word                  ([[:alnum:]]|[[:punct:]])*
string                {word}({some_space}{word})*

expon                 [Ee][-+]?{digit}+
fract                 [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fract}{expon}?)|({digit}+{expon}))

x                     {floatNum}
y                     {floatNum}
z                     {floatNum}

solid                 {space}("solid"|"SOLID"){space}
color                 {space}("color"|"COLOR"){some_space}{floatNum}{some_space}{floatNum}{some_space}{floatNum}{space}
facet                 {space}("facet"|"FACET"){space}
normal                {space}("normal"|"NORMAL"){space}
point                 {space}{x}{some_space}{y}{some_space}{z}{space}
outerloop             {space}("outer"{some_space}"loop")|("OUTER"{some_space}"LOOP"){space}
vertex                {space}("vertex"|"VERTEX"){space}
endloop               {space}("endloop"|"ENDLOOP"){space}
endfacet              {space}("endfacet"|"ENDFACET"){space}
endsolid              {space}("endsolid"|"ENDSOLID")({some_space}{word})*


 /* ------------------------------------------------------------------------ *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------ */

%option stack

%x readSolidName
%x readFacet
%x readNormal
%x readVertices
%x readVertex
%x stlError

%%

%{
    static const char* stateNames[7] =
    {
        "reading solid",
        "reading solid name",
        "reading facet",
        "reading normal",
        "reading vertices",
        "reading vertex",
        "error"
    };

    static const char* stateExpects[7] =
    {
        "'solid', 'color', 'facet' or 'endsolid'",
        "<string>",
        "'normal', 'outer loop' or 'endfacet'",
        "<float> <float> <float>",
        "'vertex' or 'endloop'",
        "<float> <float> <float>",
        ""
    };
%}


 /* ------------------------------------------------------------------------ *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------ */

 /*                      ------ Reading control header ------                */

{solid} {
    BEGIN(readSolidName);
}

<readSolidName>{string} {
    beginSolid(word::validate(YYText()));
    BEGIN(INITIAL);
}

<readSolidName>{space}\n {
    beginSolid("solid");  // Could also use solid0, solid1, ...
    ++lineNum_;
    BEGIN(INITIAL);
}

{color} {
    /* ignore 'color' */
}

{facet} {
    beginFacet();
    BEGIN(readFacet);
}

<readFacet>{normal} {
    BEGIN(readNormal);
}

<readNormal>{point} {
    /*
     skip reading normals:
     normal.x() = strtof(YYText(), &endPtr);
     normal.y() = strtof(endPtr, &endPtr);
     normal.z() = strtof(endPtr, &endPtr);
     normals_.push_back(normal);
     */
    BEGIN(readFacet);
}

<readFacet>{outerloop} {
    BEGIN(readVertices);
}

<readVertices>{vertex} {
    BEGIN(readVertex);
}

<readVertex>{space}{intNum}{space} {
    if (addVertexComponent(float(::atol(YYText()))))
    {
        BEGIN(readVertices);
    }
}

<readVertex>{space}{floatNum}{space} {
    if (addVertexComponent(::atof(YYText())))
    {
        BEGIN(readVertices);
    }
}

<readVertices>{endloop} {
    BEGIN(readFacet);
}

<readFacet>{endfacet} {
    endFacet();
    BEGIN(INITIAL);
}

{endsolid} {
}


 /* ---------------- Ignore remaining spaces and newlines ------------------ */

<*>{space} {}
<*>\n      { ++lineNum_; }


 /* ------------------- Any other characters are errors -------------------- */

<*>. {
    startError_ = YYText();
    yy_push_state(stlError);
}


 /* ---------------------------- Error handler ----------------------------- */

<stlError>.* {
    yy_pop_state();
    FatalErrorInFunction
        << "while " << stateNames[YY_START] << " on line " << lineNum_ << nl
        << "    expected " << stateExpects[YY_START]
        << " but found '" << startError_.c_str() << YYText() << "'"
        << exit(FatalError);
}


 /*  ------------------------ On EOF terminate ----------------------------  */

<<EOF>> {
    yyterminate();
}
%%


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STLReader::readAsciiFlex
(
    const fileName& filename
)
{
    IFstream is(filename);
    if (!is)
    {
        FatalErrorInFunction
            << "file " << filename << " not found"
            << exit(FatalError);
    }

    // Create with estimated number of triangles in the STL.
    // 180 bytes / triangle. For simplicity, ignore compression

    const auto fileLen = is.fileSize();

    const label nTrisEstimated =
    (
        (fileLen > 0)
      ? max(label(100), label(fileLen/180))
      : label(100)
    );


    STLAsciiParseFlex lexer(&(is.stdStream()), nTrisEstimated);
    lexer.execute();

    transfer(lexer);

    return true;
}


// ************************************************************************* //
