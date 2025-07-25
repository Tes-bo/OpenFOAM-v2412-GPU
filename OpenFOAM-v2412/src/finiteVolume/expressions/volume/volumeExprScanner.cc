
#line 1 "volumeExprScanner.rl"
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Ragel lexer interface for lemon grammar of volume field expressions

\*---------------------------------------------------------------------------*/

#include "volumeExprScanner.H"
#include "volumeExprDriver.H"
#include "volumeExprLemonParser.h"
#include "volumeExprParser.H"
#include "Enum.H"
//#include "macros.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wunused-const-variable"

// Debugging to stderr
#undef  DebugInfo
#define DebugInfo if (debug & 0x2) InfoErr

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

//- Paste token prefix
#define TOKEN_OF(T)         TOK_##T

//- An {int, c_str} enum pairing
#define TOKEN_PAIR(Name,T)  { TOKEN_OF(T), Name }

//- An {int, c_str} enum pairing for field types
#define FIELD_PAIR(Fld,T)   { TOKEN_OF(T), Fld::typeName.c_str() }


// Special handling for these known (stashed) look-back types
#define HAS_LOOKBEHIND_TOKENS
static const Enum<int> lookBehindTokenEnums
({
    TOKEN_PAIR("cellZone", CELL_ZONE), TOKEN_PAIR("cellSet", CELL_SET),
    TOKEN_PAIR("faceZone", FACE_ZONE), TOKEN_PAIR("faceSet", FACE_SET),
    #ifdef TOK_POINT_ZONE
    TOKEN_PAIR("pointZone", POINT_ZONE), TOKEN_PAIR("pointSet", POINT_SET),
    #endif
});


// Special handling of predefined method types. Eg, .x(), .y(), ...
static const Enum<int> fieldMethodEnums
({
    TOKEN_PAIR("x", CMPT_X),
    TOKEN_PAIR("y", CMPT_Y),
    TOKEN_PAIR("z", CMPT_Z),
    TOKEN_PAIR("xx", CMPT_XX),
    TOKEN_PAIR("xy", CMPT_XY),
    TOKEN_PAIR("xz", CMPT_XZ),
    TOKEN_PAIR("yx", CMPT_YX),
    TOKEN_PAIR("yy", CMPT_YY),
    TOKEN_PAIR("yz", CMPT_YZ),
    TOKEN_PAIR("zx", CMPT_ZX),
    TOKEN_PAIR("zy", CMPT_ZY),
    TOKEN_PAIR("zz", CMPT_ZZ),
    TOKEN_PAIR("ii", CMPT_II),
    TOKEN_PAIR("diag", DIAG),   /* tensors only */
    TOKEN_PAIR("T", TRANSPOSE), /* tensors only */
});


// Known field-token types
static const Enum<int>& fieldTokenEnums()
{
    static Enum<int> enums_;

    if (enums_.empty())
    {
        enums_.push_back
        ({
        #ifdef TOK_SCALAR_ID
            FIELD_PAIR(volScalarField, SCALAR_ID),
            FIELD_PAIR(volVectorField, VECTOR_ID),
            FIELD_PAIR(volTensorField, TENSOR_ID),
            FIELD_PAIR(volSymmTensorField, SYM_TENSOR_ID),
            FIELD_PAIR(volSphericalTensorField, SPH_TENSOR_ID),
        #else
            #error TOK_SCALAR_ID not defined
        #endif
        #ifdef TOK_SSCALAR_ID
            FIELD_PAIR(surfaceScalarField, SSCALAR_ID),
            FIELD_PAIR(surfaceVectorField, SVECTOR_ID),
            FIELD_PAIR(surfaceTensorField, STENSOR_ID),
            FIELD_PAIR(surfaceSymmTensorField, SSYM_TENSOR_ID),
            FIELD_PAIR(surfaceSphericalTensorField, SSPH_TENSOR_ID),
        #else
            #warning TOK_SSCALAR_ID not defined
        #endif
        #ifdef TOK_PSCALAR_ID
            FIELD_PAIR(pointScalarField, PSCALAR_ID),
            FIELD_PAIR(pointVectorField, PVECTOR_ID),
            FIELD_PAIR(pointTensorField, PTENSOR_ID),
            FIELD_PAIR(pointSymmTensorField, PSYM_TENSOR_ID),
            FIELD_PAIR(pointSphericalTensorField, PSPH_TENSOR_ID),
        #else
            #warning TOK_PSCALAR_ID not defined
        #endif
        });
    }

    return enums_;
}


// Simple compile-time function name declarations.
// Useful for handling driver-specific dispatching, or functions that
// are not universally available.
static const Enum<int> funcTokenEnums
({
#ifdef TOK_FLOOR
    TOKEN_PAIR("floor", FLOOR),
    TOKEN_PAIR("ceil", CEIL),
    TOKEN_PAIR("round", ROUND),
#endif
#ifdef TOK_HYPOT
    TOKEN_PAIR("hypot", HYPOT),
#endif

    // Already parsed as function: TOKEN_PAIR("pos", CELL_CENTRE),

    TOKEN_PAIR("point", POINT_EXPR), // Point value
    TOKEN_PAIR("face", FACE_EXPR),   // Face value, Face areaNormal

    TOKEN_PAIR("cellToPoint", CELL_TO_POINT),
    TOKEN_PAIR("cellToFace",  CELL_TO_FACE),
    TOKEN_PAIR("pointToCell", POINT_TO_CELL),

    TOKEN_PAIR("area", FACE_AREA),
    TOKEN_PAIR("vol",  CELL_VOLUME),

    TOKEN_PAIR("fpos", FACE_CENTRE),
    TOKEN_PAIR("pts", POINTS),
});

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Classifying token type based on an identifier name is indeed ugly.
//
// 1)
//   Handle special cases (eg, cellSet,...) first that have been tagged
//   as expected content with the stashed "look-behind" token.
//   Handle not-found errors here directly.
//
// 2)
//   Fallback to determining which field-type (volScalarField etc) the name
//   corresponds to.
//   Handle not-found errors by return -1.
//
static int driverTokenType
(
    const expressions::volumeExpr::parseDriver& driver_,
    const word& ident
)
{
    #ifdef HAS_LOOKBEHIND_TOKENS
    // Get stashed "look-behind" to decide what type of identifier we expect
    const int lookBehind = driver_.resetStashedTokenId();

    if (lookBehind && lookBehindTokenEnums.found(lookBehind))
    {
        bool good = false;

        switch (lookBehind)
        {
            case TOK_CELL_ZONE : good = driver_.isCellZone(ident); break;
            case TOK_CELL_SET : good = driver_.isCellSet(ident); break;

            case TOK_FACE_ZONE : good = driver_.isFaceZone(ident); break;
            case TOK_FACE_SET : good = driver_.isFaceSet(ident); break;

            #ifdef TOK_POINT_ZONE
            case TOK_POINT_ZONE : good = driver_.isPointZone(ident); break;
            case TOK_POINT_SET : good = driver_.isPointSet(ident); break;
            #endif
        }

        if (good)
        {
            return TOK_IDENTIFIER;
        }

        // Fatal
        driver_.reportFatal
        (
            "Error no " + lookBehindTokenEnums.get(lookBehind) + ": " + ident
        );

        return -2;  // Extra safety
    }
    #endif

    // Surface variables - distinguish from volume by size
    #ifdef TOK_SSCALAR_ID
    {
        const label len = driver_.mesh().nInternalFaces();

        #undef  doLocalCode
        #define doLocalCode(TokType, Type)                          \
        if (driver_.isVariable<Type>(ident, false, len))            \
        {                                                           \
            return TokType;                                         \
        }

        doLocalCode(TOK_SSCALAR_ID, scalar);
        doLocalCode(TOK_SVECTOR_ID, vector);
        doLocalCode(TOK_SSYM_TENSOR_ID, symmTensor);
        doLocalCode(TOK_SSPH_TENSOR_ID, sphericalTensor);
        doLocalCode(TOK_STENSOR_ID, tensor);
        #undef doLocalCode
    }
    #endif

    // Point variables
    #ifdef TOK_PSCALAR_ID
    {
        #undef  doLocalCode
        #define doLocalCode(TokType, Type)                          \
        if (driver_.isVariable<Type>(ident, true))                  \
        {                                                           \
            return TokType;                                         \
        }

        doLocalCode(TOK_PSCALAR_ID, scalar);
        doLocalCode(TOK_PVECTOR_ID, vector);
        doLocalCode(TOK_PSPH_TENSOR_ID, sphericalTensor);
        doLocalCode(TOK_PSYM_TENSOR_ID, symmTensor);
        doLocalCode(TOK_PTENSOR_ID, tensor);
        #undef doLocalCode
    }
    #endif

    // Volume variables
    #ifdef TOK_SCALAR_ID
    {
        #undef  doLocalCode
        #define doLocalCode(TokType, Type)                          \
        if (driver_.isVariable<Type>(ident, false))                 \
        {                                                           \
            return TokType;                                         \
        }

        doLocalCode(TOK_SCALAR_ID, scalar);
        doLocalCode(TOK_VECTOR_ID, vector);
        doLocalCode(TOK_SPH_TENSOR_ID, sphericalTensor);
        doLocalCode(TOK_SYM_TENSOR_ID, symmTensor);
        doLocalCode(TOK_TENSOR_ID, tensor);
        #undef doLocalCode
    }
    #endif

    // Check registered fields and/or disk-files
    {
        const word fieldType(driver_.getFieldClassName(ident));

        int tokType = fieldTokenEnums().lookup(fieldType, -1);

        if (tokType > 0)
        {
            return tokType;
        }
    }

    return -1;
}

} // End namespace Foam


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

#define EMIT_TOKEN(T)                                                         \
    driver_.parsePosition() = (ts-buf);                                       \
    DebugInfo<< #T << " at " << driver_.parsePosition() << nl;                \
    parser_->parse(TOKEN_OF(T));                                              \
    driver_.parsePosition() = (p-buf);


#define EMIT_VECTOR_TOKEN(X, Y, Z)                                            \
    driver_.parsePosition() = (ts-buf);                                       \
    DebugInfo<< "VECTOR at " << driver_.parsePosition() << nl;                \
    scanToken scanTok;                                                        \
    scanTok.setVector(X,Y,Z);                                                 \
    parser_->parse(TOK_VECTOR_VALUE, scanTok);                                \
    driver_.parsePosition() = (p-buf);



#line 333 "volumeExprScanner.cc"
static const int volumeExpr_start = 14;
static const int volumeExpr_first_final = 14;
static const int volumeExpr_error = 0;

static const int volumeExpr_en_main = 14;


#line 484 "volumeExprScanner.rl"



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressions::volumeExpr::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::expressions::volumeExpr::scanner::dispatch_method
(
    const parseDriver& driver_,
    word ident
) const
{
    if (ident[0] == '.')
    {
        ident.erase(0, 1);
    }

    DebugInfo
        << "Method:" << ident
        << " at " << driver_.parsePosition() << nl;

    const int methType = fieldMethodEnums.lookup(ident, -1);

    if (methType > 0)
    {
        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT);
        parser_->parse(methType);

        return true;
    }

    driver_.reportFatal("Unknown method: " + ident);
    return false;
}


bool Foam::expressions::volumeExpr::scanner::dispatch_ident
(
    const parseDriver& driver_,
    word ident
) const
{
    // Peek at stashed "look-behind". It may influence decisions
    int lookBehindTok = driver_.stashedTokenId();
    int tokType = -1;

    const bool quoted =
    (
        (ident.front() == '"' || ident.front() == '\'')
     && (ident.front() == ident.back())
    );

    if (quoted)
    {
        ident.erase(ident.size()-1);
        ident.erase(0, 1);
    }
    else
    {
        // Check for function name
        tokType = funcTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " function:"
                << parser_->tokenName(tokType) << nl;

            parser_->parse(tokType);
            return true;
        }

        #ifdef HAS_LOOKBEHIND_TOKENS
        // Specials such "cellSet" etc also reset the look-behind
        tokType = lookBehindTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " as look-behind:"
                << parser_->tokenName(tokType) << nl;

            driver_.resetStashedTokenId(tokType);
            parser_->parse(tokType);
            return true;
        }
        #endif
    }

    // Functions: scalar, vector, probably don't need others
    // - "fn:" prefix to avoid any ambiguities
    if (lookBehindTok <= 0 && ident.starts_with("fn:"))
    {
        word funcName(ident.substr(3));  // strip prefix

        do
        {
            #undef  doLocalCode
            #define doLocalCode(TokType, Type)                      \
            if (driver_.isFunction<Type>(funcName))                 \
            {                                                       \
                ident = std::move(funcName);                        \
                tokType = TokType;                                  \
                break;                                              \
            }

            #ifdef TOK_SCALAR_FUNCTION_ID
            doLocalCode(TOK_SCALAR_FUNCTION_ID, scalar);
            #endif
            #ifdef TOK_VECTOR_FUNCTION_ID
            doLocalCode(TOK_VECTOR_FUNCTION_ID, vector);
            #endif
            #undef doLocalCode
        }
        while (false);
    }

    if (tokType <= 0)
    {
        tokType = driverTokenType(driver_, ident);
    }

    if (tokType > 0)
    {
        DebugInfo
            << "Emit:" << ident << " token:"
            << parser_->tokenName(tokType) << nl;

        scanToken scanTok;
        scanTok.setWord(ident);
        parser_->parse(tokType, scanTok);

        return true;
    }


    // Not found? Attempt to strip off '.x' endings etc,
    // but not when quoted

    const auto dot = ident.rfind('.');
    const int methType =
    (
        quoted || dot == std::string::npos
      ? -1
      : fieldMethodEnums.lookup(ident.substr(dot+1), -1)
    );

    if
    (
        methType > 0
     && (tokType = driverTokenType(driver_, ident.substr(0, dot))) > 0
    )
    {
        DebugInfo
            << "Emit:" << ident.substr(0, dot).c_str() << " token:"
            << parser_->tokenName(tokType) << " with "
            << ident.substr(dot).c_str() << " token:"
            << parser_->tokenName(methType) << nl;

        // The field (before the ".")
        ident.erase(dot);

        scanToken scanTok;
        scanTok.setWord(ident);
        parser_->parse(tokType, scanTok);

        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT);
        parser_->parse(methType);

        return true;
    }

    driver_.reportFatal
    (
        "Object " + ident + " does not exist or wrong type"
    );

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::volumeExpr::scanner::process
(
    const std::string& str,
    size_t strBeg,
    size_t strLen,
    parseDriver& driver_
)
{
    // Save debug value
    const int oldDebug = debug;

    if (driver_.debugScanner()) { debug |= 0x2; }
    if (driver_.debugParser())  { debug |= 0x4; }

    if (debug & 0x6)
    {
        InfoErr
            << "Begin parse {"
            << str.substr(strBeg, strLen).c_str() << '}' << nl;
    }

    if (!parser_)
    {
        parser_ = new parser();
    }

    driver_.content(str, strBeg, strLen);

    size_t strEnd = str.length();

    if (strBeg > str.length())
    {
        strBeg = str.length();
    }
    else if (strLen != std::string::npos)
    {
        strLen += strBeg;

        if (strLen < str.length())
        {
            strEnd = strLen;
        }
    }


    parser_->start(driver_);

    // Token start/end (Ragel naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are Ragel naming
    // - buf is our own naming

    const char* buf = &(str[strBeg]);
    const char* eof = &(str[strEnd]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    
#line 599 "volumeExprScanner.cc"
	{
	cs = volumeExpr_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 740 "volumeExprScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 611 "volumeExprScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 358 "volumeExprScanner.rl"
	{te = p+1;{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr4:
#line 358 "volumeExprScanner.rl"
	{te = p+1;{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr5:
#line 333 "volumeExprScanner.rl"
	{{p = ((te))-1;}{
        // Emit number
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        scanToken scanTok;
        scanTok.setScalar(0);
        if (readScalar(std::string(ts, te-ts), scanTok.scalarValue))
        {
            parser_->parse(TOKEN_OF(NUMBER), scanTok);
        }
        else
        {
            driver_.reportFatal
            (
                "Error parsing number: " + std::string(ts, te-ts)
            );
        }

        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr8:
#line 406 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(EQUAL); }}
	goto st14;
tr9:
#line 358 "volumeExprScanner.rl"
	{{p = ((te))-1;}{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr11:
#line 460 "volumeExprScanner.rl"
	{{p = ((te))-1;}{ EMIT_TOKEN(TENSOR); }}
	goto st14;
tr13:
#line 471 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(IDENTITY_TENSOR); }}
	goto st14;
tr14:
#line 459 "volumeExprScanner.rl"
	{{p = ((te))-1;}{ EMIT_TOKEN(VECTOR); }}
	goto st14;
tr16:
#line 468 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_VECTOR_TOKEN(1,0,0); }}
	goto st14;
tr17:
#line 469 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_VECTOR_TOKEN(0,1,0); }}
	goto st14;
tr18:
#line 470 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_VECTOR_TOKEN(0,0,1); }}
	goto st14;
tr19:
#line 409 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LOR); }}
	goto st14;
tr23:
#line 391 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PERCENT); }}
	goto st14;
tr26:
#line 392 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LPAREN); }}
	goto st14;
tr27:
#line 393 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(RPAREN); }}
	goto st14;
tr28:
#line 394 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(TIMES); }}
	goto st14;
tr29:
#line 395 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PLUS); }}
	goto st14;
tr30:
#line 397 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COMMA); }}
	goto st14;
tr31:
#line 396 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(MINUS); }}
	goto st14;
tr33:
#line 399 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(DIVIDE); }}
	goto st14;
tr35:
#line 401 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COLON); }}
	goto st14;
tr39:
#line 400 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(QUESTION); }}
	goto st14;
tr41:
#line 412 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(BIT_XOR); }}
	goto st14;
tr58:
#line 385 "volumeExprScanner.rl"
	{te = p;p--;}
	goto st14;
tr59:
#line 390 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LNOT); }}
	goto st14;
tr60:
#line 407 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(NOT_EQUAL); }}
	goto st14;
tr61:
#line 410 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(BIT_AND); }}
	goto st14;
tr62:
#line 408 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LAND); }}
	goto st14;
tr63:
#line 398 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(DOT); }}
	goto st14;
tr66:
#line 333 "volumeExprScanner.rl"
	{te = p;p--;{
        // Emit number
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        scanToken scanTok;
        scanTok.setScalar(0);
        if (readScalar(std::string(ts, te-ts), scanTok.scalarValue))
        {
            parser_->parse(TOKEN_OF(NUMBER), scanTok);
        }
        else
        {
            driver_.reportFatal
            (
                "Error parsing number: " + std::string(ts, te-ts)
            );
        }

        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr68:
#line 365 "volumeExprScanner.rl"
	{te = p;p--;{
        // Tokenized ".method" - dispatch '.' and "method" separately
        driver_.parsePosition() = (ts-buf);
        dispatch_method(driver_, word(ts+1, te-ts-1, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr69:
#line 402 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LESS); }}
	goto st14;
tr70:
#line 403 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LESS_EQUAL); }}
	goto st14;
tr71:
#line 404 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(GREATER); }}
	goto st14;
tr72:
#line 405 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(GREATER_EQUAL); }}
	goto st14;
tr73:
#line 358 "volumeExprScanner.rl"
	{te = p;p--;{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st14;
tr75:
#line 1 "NONE"
	{	switch( act ) {
	case 26:
	{{p = ((te))-1;} EMIT_TOKEN(PI); }
	break;
	case 27:
	{{p = ((te))-1;} EMIT_TOKEN(DEG_TO_RAD); }
	break;
	case 28:
	{{p = ((te))-1;} EMIT_TOKEN(RAD_TO_DEG); }
	break;
	case 29:
	{{p = ((te))-1;} EMIT_TOKEN(EXP); }
	break;
	case 31:
	{{p = ((te))-1;} EMIT_TOKEN(LOG10); }
	break;
	case 32:
	{{p = ((te))-1;} EMIT_TOKEN(POW); }
	break;
	case 34:
	{{p = ((te))-1;} EMIT_TOKEN(SQRT); }
	break;
	case 35:
	{{p = ((te))-1;} EMIT_TOKEN(CBRT); }
	break;
	case 39:
	{{p = ((te))-1;} EMIT_TOKEN(ASIN); }
	break;
	case 40:
	{{p = ((te))-1;} EMIT_TOKEN(ACOS); }
	break;
	case 42:
	{{p = ((te))-1;} EMIT_TOKEN(ATAN2); }
	break;
	case 43:
	{{p = ((te))-1;} EMIT_TOKEN(SINH); }
	break;
	case 44:
	{{p = ((te))-1;} EMIT_TOKEN(COSH); }
	break;
	case 45:
	{{p = ((te))-1;} EMIT_TOKEN(TANH); }
	break;
	case 47:
	{{p = ((te))-1;} EMIT_TOKEN(MAGSQR); }
	break;
	case 50:
	{{p = ((te))-1;} EMIT_TOKEN(POS0); }
	break;
	case 51:
	{{p = ((te))-1;} EMIT_TOKEN(NEG0); }
	break;
	case 52:
	{{p = ((te))-1;} EMIT_TOKEN(SIGN); }
	break;
	case 53:
	{{p = ((te))-1;} EMIT_TOKEN(MIN); }
	break;
	case 54:
	{{p = ((te))-1;} EMIT_TOKEN(MAX); }
	break;
	case 55:
	{{p = ((te))-1;} EMIT_TOKEN(AVERAGE); }
	break;
	case 56:
	{{p = ((te))-1;} EMIT_TOKEN(SUM); }
	break;
	case 57:
	{{p = ((te))-1;} EMIT_TOKEN(WEIGHT_AVERAGE); }
	break;
	case 58:
	{{p = ((te))-1;} EMIT_TOKEN(WEIGHT_SUM); }
	break;
	case 59:
	{{p = ((te))-1;} EMIT_TOKEN(RAND); }
	break;
	case 60:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL); }
	break;
	case 63:
	{{p = ((te))-1;} EMIT_TOKEN(SYM_TENSOR); }
	break;
	case 64:
	{{p = ((te))-1;} EMIT_TOKEN(SPH_TENSOR); }
	break;
	case 65:
	{{p = ((te))-1;} EMIT_TOKEN(LTRUE); }
	break;
	case 66:
	{{p = ((te))-1;} EMIT_TOKEN(LFALSE); }
	break;
	case 67:
	{{p = ((te))-1;} EMIT_TOKEN(ZERO); }
	break;
	case 72:
	{{p = ((te))-1;} EMIT_TOKEN(ARG); }
	break;
	case 73:
	{{p = ((te))-1;} EMIT_TOKEN(TIME); }
	break;
	case 74:
	{{p = ((te))-1;} EMIT_TOKEN(DELTA_T); }
	break;
	case 75:
	{{p = ((te))-1;}
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }
	break;
	}
	}
	goto st14;
tr91:
#line 434 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(ATAN); }}
	goto st14;
tr106:
#line 430 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(COS); }}
	goto st14;
tr129:
#line 423 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LOG); }}
	goto st14;
tr136:
#line 439 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(MAG); }}
	goto st14;
tr143:
#line 443 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NEG); }}
	goto st14;
tr149:
#line 442 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(POS); }}
	goto st14;
tr168:
#line 429 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SIN); }}
	goto st14;
tr184:
#line 426 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SQR); }}
	goto st14;
tr200:
#line 431 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TAN); }}
	goto st14;
tr206:
#line 460 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TENSOR); }}
	goto st14;
tr217:
#line 459 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(VECTOR); }}
	goto st14;
st14:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 1 "NONE"
	{ts = p;}
#line 1000 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 32: goto st15;
		case 33: goto st16;
		case 34: goto st1;
		case 37: goto tr23;
		case 38: goto st17;
		case 39: goto st3;
		case 40: goto tr26;
		case 41: goto tr27;
		case 42: goto tr28;
		case 43: goto tr29;
		case 44: goto tr30;
		case 45: goto tr31;
		case 46: goto st18;
		case 47: goto tr33;
		case 58: goto tr35;
		case 60: goto st23;
		case 61: goto st7;
		case 62: goto st24;
		case 63: goto tr39;
		case 90: goto st27;
		case 94: goto tr41;
		case 95: goto st25;
		case 97: goto st30;
		case 98: goto st44;
		case 99: goto st47;
		case 100: goto st52;
		case 101: goto st62;
		case 102: goto st64;
		case 108: goto st69;
		case 109: goto st73;
		case 110: goto st79;
		case 112: goto st82;
		case 114: goto st85;
		case 115: goto st93;
		case 116: goto st121;
		case 118: goto st133;
		case 119: goto st139;
		case 124: goto st13;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st15;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 89 ) {
			if ( 103 <= (*p) && (*p) <= 122 )
				goto st25;
		} else if ( (*p) >= 65 )
			goto st25;
	} else
		goto tr34;
	goto st0;
st0:
cs = 0;
	goto _out;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) == 32 )
		goto st15;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st15;
	goto tr58;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	if ( (*p) == 61 )
		goto tr60;
	goto tr59;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 34 )
		goto st0;
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( (*p) == 34 )
		goto tr2;
	goto st2;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) == 38 )
		goto tr62;
	goto tr61;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 39 )
		goto st0;
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 39 )
		goto tr4;
	goto st4;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr64;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st21;
	} else
		goto st21;
	goto tr63;
tr64:
#line 1 "NONE"
	{te = p+1;}
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 1128 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr64;
	goto tr66;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 43: goto st6;
		case 45: goto st6;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st20;
	goto tr5;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st20;
	goto tr5;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st20;
	goto tr66;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st21;
	} else if ( (*p) >= 65 )
		goto st21;
	goto tr68;
tr34:
#line 1 "NONE"
	{te = p+1;}
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
#line 1179 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr64;
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr34;
	goto tr66;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 61 )
		goto tr70;
	goto tr69;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 61 )
		goto tr8;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 61 )
		goto tr72;
	goto tr71;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
tr74:
#line 1 "NONE"
	{te = p+1;}
#line 358 "volumeExprScanner.rl"
	{act = 75;}
	goto st26;
tr78:
#line 1 "NONE"
	{te = p+1;}
#line 467 "volumeExprScanner.rl"
	{act = 67;}
	goto st26;
tr85:
#line 1 "NONE"
	{te = p+1;}
#line 433 "volumeExprScanner.rl"
	{act = 40;}
	goto st26;
tr86:
#line 1 "NONE"
	{te = p+1;}
#line 472 "volumeExprScanner.rl"
	{act = 72;}
	goto st26;
tr88:
#line 1 "NONE"
	{te = p+1;}
#line 432 "volumeExprScanner.rl"
	{act = 39;}
	goto st26;
tr92:
#line 1 "NONE"
	{te = p+1;}
#line 435 "volumeExprScanner.rl"
	{act = 42;}
	goto st26;
tr97:
#line 1 "NONE"
	{te = p+1;}
#line 451 "volumeExprScanner.rl"
	{act = 55;}
	goto st26;
tr100:
#line 1 "NONE"
	{te = p+1;}
#line 458 "volumeExprScanner.rl"
	{act = 60;}
	goto st26;
tr104:
#line 1 "NONE"
	{te = p+1;}
#line 428 "volumeExprScanner.rl"
	{act = 35;}
	goto st26;
tr107:
#line 1 "NONE"
	{te = p+1;}
#line 437 "volumeExprScanner.rl"
	{act = 44;}
	goto st26;
tr115:
#line 1 "NONE"
	{te = p+1;}
#line 420 "volumeExprScanner.rl"
	{act = 27;}
	goto st26;
tr118:
#line 1 "NONE"
	{te = p+1;}
#line 474 "volumeExprScanner.rl"
	{act = 74;}
	goto st26;
tr120:
#line 1 "NONE"
	{te = p+1;}
#line 422 "volumeExprScanner.rl"
	{act = 29;}
	goto st26;
tr125:
#line 1 "NONE"
	{te = p+1;}
#line 466 "volumeExprScanner.rl"
	{act = 66;}
	goto st26;
tr131:
#line 1 "NONE"
	{te = p+1;}
#line 424 "volumeExprScanner.rl"
	{act = 31;}
	goto st26;
tr135:
#line 1 "NONE"
	{te = p+1;}
#line 450 "volumeExprScanner.rl"
	{act = 54;}
	goto st26;
tr139:
#line 1 "NONE"
	{te = p+1;}
#line 440 "volumeExprScanner.rl"
	{act = 47;}
	goto st26;
tr140:
#line 1 "NONE"
	{te = p+1;}
#line 449 "volumeExprScanner.rl"
	{act = 53;}
	goto st26;
tr144:
#line 1 "NONE"
	{te = p+1;}
#line 445 "volumeExprScanner.rl"
	{act = 51;}
	goto st26;
tr145:
#line 1 "NONE"
	{te = p+1;}
#line 419 "volumeExprScanner.rl"
	{act = 26;}
	goto st26;
tr148:
#line 1 "NONE"
	{te = p+1;}
#line 425 "volumeExprScanner.rl"
	{act = 32;}
	goto st26;
tr150:
#line 1 "NONE"
	{te = p+1;}
#line 444 "volumeExprScanner.rl"
	{act = 50;}
	goto st26;
tr158:
#line 1 "NONE"
	{te = p+1;}
#line 421 "volumeExprScanner.rl"
	{act = 28;}
	goto st26;
tr159:
#line 1 "NONE"
	{te = p+1;}
#line 455 "volumeExprScanner.rl"
	{act = 59;}
	goto st26;
tr167:
#line 1 "NONE"
	{te = p+1;}
#line 446 "volumeExprScanner.rl"
	{act = 52;}
	goto st26;
tr169:
#line 1 "NONE"
	{te = p+1;}
#line 436 "volumeExprScanner.rl"
	{act = 43;}
	goto st26;
tr182:
#line 1 "NONE"
	{te = p+1;}
#line 462 "volumeExprScanner.rl"
	{act = 64;}
	goto st26;
tr185:
#line 1 "NONE"
	{te = p+1;}
#line 427 "volumeExprScanner.rl"
	{act = 34;}
	goto st26;
tr186:
#line 1 "NONE"
	{te = p+1;}
#line 452 "volumeExprScanner.rl"
	{act = 56;}
	goto st26;
tr194:
#line 1 "NONE"
	{te = p+1;}
#line 461 "volumeExprScanner.rl"
	{act = 63;}
	goto st26;
tr201:
#line 1 "NONE"
	{te = p+1;}
#line 438 "volumeExprScanner.rl"
	{act = 45;}
	goto st26;
tr209:
#line 1 "NONE"
	{te = p+1;}
#line 473 "volumeExprScanner.rl"
	{act = 73;}
	goto st26;
tr211:
#line 1 "NONE"
	{te = p+1;}
#line 465 "volumeExprScanner.rl"
	{act = 65;}
	goto st26;
tr231:
#line 1 "NONE"
	{te = p+1;}
#line 453 "volumeExprScanner.rl"
	{act = 57;}
	goto st26;
tr233:
#line 1 "NONE"
	{te = p+1;}
#line 454 "volumeExprScanner.rl"
	{act = 58;}
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
#line 1440 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr75;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st28;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto st29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto tr78;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 99: goto st31;
		case 114: goto st33;
		case 115: goto st34;
		case 116: goto st36;
		case 118: goto st39;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto tr86;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 105: goto st35;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto tr88;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto st38;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 46: goto tr74;
		case 50: goto tr92;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr91;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st40;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto st41;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st42;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st43;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st45;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st46;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 108: goto tr100;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 98: goto st48;
		case 111: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto st49;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 116: goto tr104;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto st51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 104: goto tr107;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr106;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st54;
		case 108: goto st59;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 46: goto tr74;
		case 84: goto st55;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st56;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 46: goto tr74;
		case 82: goto st57;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st58;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 100: goto tr115;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 116: goto st60;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st61;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 46: goto tr74;
		case 84: goto tr118;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 120: goto st63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 112: goto tr120;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st65;
		case 110: goto tr122;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 108: goto st66;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto tr125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
tr122:
#line 1 "NONE"
	{te = p+1;}
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
#line 2207 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr74;
		case 58: goto st8;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 95 )
		goto st25;
	if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st25;
	} else if ( (*p) >= 65 )
		goto st25;
	goto tr9;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st70;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 46: goto tr74;
		case 49: goto st72;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr129;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 46: goto tr74;
		case 48: goto tr131;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st74;
		case 105: goto st78;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st75;
		case 120: goto tr135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 46: goto tr74;
		case 83: goto st76;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr136;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 113: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto tr139;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto tr140;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st80;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	switch( (*p) ) {
		case 46: goto tr74;
		case 48: goto tr144;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr143;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 105: goto tr145;
		case 111: goto st83;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto st84;
		case 119: goto tr148;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 46: goto tr74;
		case 48: goto tr150;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr149;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st86;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 100: goto st87;
		case 110: goto st92;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	switch( (*p) ) {
		case 46: goto tr74;
		case 84: goto st88;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st89;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 46: goto tr74;
		case 68: goto st90;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st91;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto tr158;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 100: goto tr159;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 105: goto st94;
		case 112: goto st97;
		case 113: goto st110;
		case 117: goto st112;
		case 121: goto st113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st95;
		case 110: goto st96;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto tr167;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 104: goto tr169;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr168;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 104: goto st98;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st99;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto st100;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 105: goto st101;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 99: goto st102;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st103;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 108: goto st104;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
	switch( (*p) ) {
		case 46: goto tr74;
		case 84: goto st105;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st106;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto st107;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto st108;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st109;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto tr182;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto st111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 116: goto tr185;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr184;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 109: goto tr186;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 109: goto st114;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 109: goto st115;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	switch( (*p) ) {
		case 46: goto tr74;
		case 84: goto st116;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st117;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto st119;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st120;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto tr194;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st122;
		case 101: goto st124;
		case 105: goto st129;
		case 114: goto st131;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto st123;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 104: goto tr201;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr200;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 110: goto st125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 115: goto st126;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st127;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto tr205;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
tr205:
#line 1 "NONE"
	{te = p+1;}
	goto st128;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
#line 3317 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr74;
		case 58: goto st9;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr206;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 58 )
		goto st10;
	goto tr11;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 73 )
		goto tr13;
	goto tr11;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 109: goto st130;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto tr209;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 117: goto st132;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto tr211;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st134;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 99: goto st135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 116: goto st136;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 111: goto st137;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st137:
	if ( ++p == pe )
		goto _test_eof137;
case 137:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto tr216;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
tr216:
#line 1 "NONE"
	{te = p+1;}
	goto st138;
st138:
	if ( ++p == pe )
		goto _test_eof138;
case 138:
#line 3516 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr74;
		case 58: goto st11;
		case 95: goto tr74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr217;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 58 )
		goto st12;
	goto tr14;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 120: goto tr16;
		case 121: goto tr17;
		case 122: goto tr18;
	}
	goto tr14;
st139:
	if ( ++p == pe )
		goto _test_eof139;
case 139:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st140;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st140:
	if ( ++p == pe )
		goto _test_eof140;
case 140:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 105: goto st141;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st141:
	if ( ++p == pe )
		goto _test_eof141;
case 141:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st142;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st142:
	if ( ++p == pe )
		goto _test_eof142;
case 142:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 104: goto st143;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st143:
	if ( ++p == pe )
		goto _test_eof143;
case 143:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 116: goto st144;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st144:
	if ( ++p == pe )
		goto _test_eof144;
case 144:
	switch( (*p) ) {
		case 46: goto tr74;
		case 65: goto st145;
		case 83: goto st151;
		case 95: goto tr74;
	}
	if ( (*p) < 66 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st145:
	if ( ++p == pe )
		goto _test_eof145;
case 145:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 118: goto st146;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st146:
	if ( ++p == pe )
		goto _test_eof146;
case 146:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto st147;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st147:
	if ( ++p == pe )
		goto _test_eof147;
case 147:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 114: goto st148;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st148:
	if ( ++p == pe )
		goto _test_eof148;
case 148:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 97: goto st149;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st149:
	if ( ++p == pe )
		goto _test_eof149;
case 149:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 103: goto st150;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st150:
	if ( ++p == pe )
		goto _test_eof150;
case 150:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 101: goto tr231;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st151:
	if ( ++p == pe )
		goto _test_eof151;
case 151:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 117: goto st152;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st152:
	if ( ++p == pe )
		goto _test_eof152;
case 152:
	switch( (*p) ) {
		case 46: goto tr74;
		case 95: goto tr74;
		case 109: goto tr233;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr74;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr74;
	} else
		goto tr74;
	goto tr73;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 124 )
		goto tr19;
	goto st0;
	}
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof148: cs = 148; goto _test_eof; 
	_test_eof149: cs = 149; goto _test_eof; 
	_test_eof150: cs = 150; goto _test_eof; 
	_test_eof151: cs = 151; goto _test_eof; 
	_test_eof152: cs = 152; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 15: goto tr58;
	case 16: goto tr59;
	case 17: goto tr61;
	case 18: goto tr63;
	case 19: goto tr66;
	case 5: goto tr5;
	case 6: goto tr5;
	case 20: goto tr66;
	case 21: goto tr68;
	case 22: goto tr66;
	case 23: goto tr69;
	case 24: goto tr71;
	case 25: goto tr73;
	case 26: goto tr75;
	case 27: goto tr73;
	case 28: goto tr73;
	case 29: goto tr73;
	case 30: goto tr73;
	case 31: goto tr73;
	case 32: goto tr73;
	case 33: goto tr73;
	case 34: goto tr73;
	case 35: goto tr73;
	case 36: goto tr73;
	case 37: goto tr73;
	case 38: goto tr91;
	case 39: goto tr73;
	case 40: goto tr73;
	case 41: goto tr73;
	case 42: goto tr73;
	case 43: goto tr73;
	case 44: goto tr73;
	case 45: goto tr73;
	case 46: goto tr73;
	case 47: goto tr73;
	case 48: goto tr73;
	case 49: goto tr73;
	case 50: goto tr73;
	case 51: goto tr106;
	case 52: goto tr73;
	case 53: goto tr73;
	case 54: goto tr73;
	case 55: goto tr73;
	case 56: goto tr73;
	case 57: goto tr73;
	case 58: goto tr73;
	case 59: goto tr73;
	case 60: goto tr73;
	case 61: goto tr73;
	case 62: goto tr73;
	case 63: goto tr73;
	case 64: goto tr73;
	case 65: goto tr73;
	case 66: goto tr73;
	case 67: goto tr73;
	case 68: goto tr73;
	case 8: goto tr9;
	case 69: goto tr73;
	case 70: goto tr73;
	case 71: goto tr129;
	case 72: goto tr73;
	case 73: goto tr73;
	case 74: goto tr73;
	case 75: goto tr136;
	case 76: goto tr73;
	case 77: goto tr73;
	case 78: goto tr73;
	case 79: goto tr73;
	case 80: goto tr73;
	case 81: goto tr143;
	case 82: goto tr73;
	case 83: goto tr73;
	case 84: goto tr149;
	case 85: goto tr73;
	case 86: goto tr73;
	case 87: goto tr73;
	case 88: goto tr73;
	case 89: goto tr73;
	case 90: goto tr73;
	case 91: goto tr73;
	case 92: goto tr73;
	case 93: goto tr73;
	case 94: goto tr73;
	case 95: goto tr73;
	case 96: goto tr168;
	case 97: goto tr73;
	case 98: goto tr73;
	case 99: goto tr73;
	case 100: goto tr73;
	case 101: goto tr73;
	case 102: goto tr73;
	case 103: goto tr73;
	case 104: goto tr73;
	case 105: goto tr73;
	case 106: goto tr73;
	case 107: goto tr73;
	case 108: goto tr73;
	case 109: goto tr73;
	case 110: goto tr73;
	case 111: goto tr184;
	case 112: goto tr73;
	case 113: goto tr73;
	case 114: goto tr73;
	case 115: goto tr73;
	case 116: goto tr73;
	case 117: goto tr73;
	case 118: goto tr73;
	case 119: goto tr73;
	case 120: goto tr73;
	case 121: goto tr73;
	case 122: goto tr73;
	case 123: goto tr200;
	case 124: goto tr73;
	case 125: goto tr73;
	case 126: goto tr73;
	case 127: goto tr73;
	case 128: goto tr206;
	case 9: goto tr11;
	case 10: goto tr11;
	case 129: goto tr73;
	case 130: goto tr73;
	case 131: goto tr73;
	case 132: goto tr73;
	case 133: goto tr73;
	case 134: goto tr73;
	case 135: goto tr73;
	case 136: goto tr73;
	case 137: goto tr73;
	case 138: goto tr217;
	case 11: goto tr14;
	case 12: goto tr14;
	case 139: goto tr73;
	case 140: goto tr73;
	case 141: goto tr73;
	case 142: goto tr73;
	case 143: goto tr73;
	case 144: goto tr73;
	case 145: goto tr73;
	case 146: goto tr73;
	case 147: goto tr73;
	case 148: goto tr73;
	case 149: goto tr73;
	case 150: goto tr73;
	case 151: goto tr73;
	case 152: goto tr73;
	}
	}

	_out: {}
	}

#line 742 "volumeExprScanner.rl"
   /* ^^^ FSM execution here ^^^ */;

    if (0 == cs)
    {
        driver_.reportFatal("Parse error while scanning", (p-buf));
    }

    if (p != eof)
    {
        driver_.reportFatal("Parsing failed with remaining content", (p-buf));
    }

    // Terminate parser execution
    parser_->parse(0);
    parser_->stop();

    if (debug & 0x6)
    {
        InfoErr<< "Done parse." << nl;
    }

    // Restore debug value
    debug = oldDebug;

    return true;
}


// ************************************************************************* //
