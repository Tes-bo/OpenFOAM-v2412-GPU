/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Ragel lexer interface for lemon grammar for patch expressions

\*---------------------------------------------------------------------------*/

#include "exprScanToken.H"
#include "fieldExprScanner.H"
#include "fieldExprDriver.H"
#include "fieldExprLemonParser.h"
#include "fieldExprParser.H"
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

// No known look-back types
#undef HAS_LOOKBEHIND_TOKENS

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
    const expressions::fieldExpr::parseDriver& driver_,
    const word& ident
)
{
    #ifdef HAS_LOOKBEHIND_TOKENS
    // Get stashed "look-behind" to decide what type of identifier we expect
    #endif

    // Field variables
    #ifdef TOK_SCALAR_ID
    {
        #undef  doLocalCode
        #define doLocalCode(TokType, Type)                          \
        if (driver_.isLocalVariable<Type>(ident, false))            \
        {                                                           \
            return TokType;                                         \
        }

        doLocalCode(TOK_SCALAR_ID, scalar);
        doLocalCode(TOK_VECTOR_ID, vector);
        doLocalCode(TOK_SYM_TENSOR_ID, symmTensor);
        doLocalCode(TOK_SPH_TENSOR_ID, sphericalTensor);
        doLocalCode(TOK_TENSOR_ID, tensor);
        // Not tested: doLocalCode(TOK_BOOL_ID, bool);
        #undef doLocalCode
    }
    #endif

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


%%{
    machine fieldExpr;
    write   data;

    action emit_number {
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
    }

    action emit_ident {
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }

    action emit_method {
        // Tokenized ".method" - dispatch '.' and "method" separately
        driver_.parsePosition() = (ts-buf);
        dispatch_method(driver_, word(ts+1, te-ts-1, false));
        driver_.parsePosition() = (p-buf);
    }

    decimal = ((digit* '.' digit+) | (digit+ '.'?)) ;
    number  = ((digit+ | decimal) ([Ee][\-+]? digit+)?) ;
    identifier = ((alpha|'_') . ((alnum|[._])**)) ;
    dquoted = '"' [^\"]+ '"' ;
    squoted = "'" [^\']+ "'" ;

    ## Allow 'fn:' prefix for function identifier
    ident = ('fn:')? identifier ;

    ## ===========
    ## The scanner
    ## ===========
    main := |*
        space*;

    number => emit_number;

    ## Operators
    '!'  =>{ EMIT_TOKEN(LNOT); };
    '%'  =>{ EMIT_TOKEN(PERCENT); };
    '('  =>{ EMIT_TOKEN(LPAREN); };
    ')'  =>{ EMIT_TOKEN(RPAREN); };
    '*'  =>{ EMIT_TOKEN(TIMES); };
    '+'  =>{ EMIT_TOKEN(PLUS); };
    '-'  =>{ EMIT_TOKEN(MINUS); };
    ','  =>{ EMIT_TOKEN(COMMA); };
    '.'  =>{ EMIT_TOKEN(DOT); };
    '/'  =>{ EMIT_TOKEN(DIVIDE); };
    '?'  =>{ EMIT_TOKEN(QUESTION); };
    ':'  =>{ EMIT_TOKEN(COLON); };
    '<'  =>{ EMIT_TOKEN(LESS); };
    '<=' =>{ EMIT_TOKEN(LESS_EQUAL); };
    '>'  =>{ EMIT_TOKEN(GREATER); };
    '>=' =>{ EMIT_TOKEN(GREATER_EQUAL); };
    '==' =>{ EMIT_TOKEN(EQUAL); };
    '!=' =>{ EMIT_TOKEN(NOT_EQUAL); };
    '&&' =>{ EMIT_TOKEN(LAND); };
    '||' =>{ EMIT_TOKEN(LOR); };
    '&'  =>{ EMIT_TOKEN(BIT_AND); };
## Not needed?  '|'  =>{ EMIT_TOKEN(BIT_OR); };
    '^'  =>{ EMIT_TOKEN(BIT_XOR); };

    ## Some '.method' - Error if unknown
    '.' alpha+ => emit_method;


    ## Regular functions
    "pi"        =>{ EMIT_TOKEN(PI); };
    "degToRad"  =>{ EMIT_TOKEN(DEG_TO_RAD); };
    "radToDeg"  =>{ EMIT_TOKEN(RAD_TO_DEG); };
    "exp"       =>{ EMIT_TOKEN(EXP); };
    "log"       =>{ EMIT_TOKEN(LOG); };
    "log10"     =>{ EMIT_TOKEN(LOG10); };
    "pow"       =>{ EMIT_TOKEN(POW); };
    "sqr"       =>{ EMIT_TOKEN(SQR); };
    "sqrt"      =>{ EMIT_TOKEN(SQRT); };
    "cbrt"      =>{ EMIT_TOKEN(CBRT); };
    "sin"       =>{ EMIT_TOKEN(SIN); };
    "cos"       =>{ EMIT_TOKEN(COS); };
    "tan"       =>{ EMIT_TOKEN(TAN); };
    "asin"      =>{ EMIT_TOKEN(ASIN); };
    "acos"      =>{ EMIT_TOKEN(ACOS); };
    "atan"      =>{ EMIT_TOKEN(ATAN); };
    "atan2"     =>{ EMIT_TOKEN(ATAN2); };
    "sinh"      =>{ EMIT_TOKEN(SINH); };
    "cosh"      =>{ EMIT_TOKEN(COSH); };
    "tanh"      =>{ EMIT_TOKEN(TANH); };
    "mag"       =>{ EMIT_TOKEN(MAG); };
    "magSqr"    =>{ EMIT_TOKEN(MAGSQR); };

    "pos"       =>{ EMIT_TOKEN(POS); };
    "neg"       =>{ EMIT_TOKEN(NEG); };
    "pos0"      =>{ EMIT_TOKEN(POS0); };
    "neg0"      =>{ EMIT_TOKEN(NEG0); };
    "sign"      =>{ EMIT_TOKEN(SIGN); };

    ## Reductions, or other special functions
    "min"       =>{ EMIT_TOKEN(MIN); };
    "max"       =>{ EMIT_TOKEN(MAX); };
    "average"   =>{ EMIT_TOKEN(AVERAGE); };
    "sum"       =>{ EMIT_TOKEN(SUM); };
    "rand"      =>{ EMIT_TOKEN(RAND); };

    ## Types
    "bool"      =>{ EMIT_TOKEN(BOOL); };
    "vector"    =>{ EMIT_TOKEN(VECTOR); };
    "tensor"    =>{ EMIT_TOKEN(TENSOR); };
    "symmTensor" =>{ EMIT_TOKEN(SYM_TENSOR); };
    "sphericalTensor" =>{ EMIT_TOKEN(SPH_TENSOR); };

    ## Single value (constants, etc)
    "true"      =>{ EMIT_TOKEN(LTRUE); };
    "false"     =>{ EMIT_TOKEN(LFALSE); };
    "Zero"      =>{ EMIT_TOKEN(ZERO); };
    "vector::x" =>{ EMIT_VECTOR_TOKEN(1,0,0); };
    "vector::y" =>{ EMIT_VECTOR_TOKEN(0,1,0); };
    "vector::z" =>{ EMIT_VECTOR_TOKEN(0,0,1); };
    "tensor::I" =>{ EMIT_TOKEN(IDENTITY_TENSOR); };
    "arg"       =>{ EMIT_TOKEN(ARG); };
    "time"      =>{ EMIT_TOKEN(TIME); };
    "deltaT"    =>{ EMIT_TOKEN(DELTA_T); };

    ## Identifier (field, etc - error if unknown)
    ## Handle 'bare' names and single/double quoted ones
    ident       => emit_ident;
    dquoted     => emit_ident;
    squoted     => emit_ident;

    space*;
    *|;
}%%


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressions::fieldExpr::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::expressions::fieldExpr::scanner::dispatch_method
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


bool Foam::expressions::fieldExpr::scanner::dispatch_ident
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

bool Foam::expressions::fieldExpr::scanner::process
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
    %%{write init;}%%   /* ^^^ FSM initialization here ^^^ */;

    %%{write exec;}%%   /* ^^^ FSM execution here ^^^ */;

    if (%%{write error;}%% == cs)
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
