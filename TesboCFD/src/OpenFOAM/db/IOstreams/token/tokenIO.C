/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "token.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
template<class OS>
static OS& printTokenInfo(OS& os, const token& tok)
{
    os  << "on line " << tok.lineNumber() << ": ";

    switch (tok.type())
    {
        case token::tokenType::UNDEFINED:
            os  << "undefined token";
        break;

        case token::tokenType::BOOL:
            os  << "bool '" << (tok.boolToken() ? "true" : "false") << '\'';
        break;

        case token::tokenType::FLAG:
            os  << "flag '" << int(tok.flagToken()) << '\'';
        break;

        case token::tokenType::PUNCTUATION:
            os  << "punctuation '" << tok.pToken() << '\'';
        break;

        // case token::tokenType::INT32:
        //     os  << "int32 " << tok.int32Token();
        // break;
        //
        // case token::tokenType::INT64:
        //     os  << "int64 " << tok.int64Token();
        // break;

        case token::tokenType::LABEL:
            os  << "label " << tok.labelToken();
        break;

        case token::tokenType::FLOAT:
            os  << "float " << tok.floatToken();
        break;

        case token::tokenType::DOUBLE:
            os  << "double " << tok.doubleToken();
        break;

        case token::tokenType::WORD:
            os  << "word '" << tok.wordToken() << '\'';
        break;

        case token::tokenType::DIRECTIVE:
            os  << "directive '" << tok.wordToken() << '\'';
        break;

        case token::tokenType::STRING:
            os  << "string " << tok.stringToken();
        break;

        case token::tokenType::EXPRESSION:
            os  << "expression " << tok.stringToken();
        break;

        case token::tokenType::VARIABLE:
            os  << "variable " << tok.stringToken();
        break;

        case token::tokenType::VERBATIM:
            os  << "verbatim " << tok.stringToken();
        break;

        case token::tokenType::CHAR_DATA:
            os  << "char_data " << tok.stringToken();
        break;

        case token::tokenType::COMPOUND:
        {
            if (tok.compoundToken().pending())
            {
                os  << "pending ";
            }
            if (tok.compoundToken().moved())
            {
                os  << "moved ";
            }
            os  << "compound of type "
                << tok.compoundToken().type();
        }
        break;

        case token::tokenType::ERROR:
            os  << "error";
        break;

        default:
            os  << "unknown token type '" << int(tok.type()) << '\'';
            break;
    }

    return os;
}
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::token::token(Istream& is)
:
    token()
{
    is.read(*this);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::token::name(const token::tokenType tokType)
{
    switch (tokType)
    {
        case token::tokenType::UNDEFINED: return "undefined";
        case token::tokenType::BOOL: return "bool";
        case token::tokenType::FLAG: return "flag";
        case token::tokenType::PUNCTUATION: return "punctuation";

        // case token::tokenType::INT32: return "int32";
        // case token::tokenType::INT64: return "int64";
        case token::tokenType::LABEL: return "label";
        case token::tokenType::FLOAT: return "float";
        case token::tokenType::DOUBLE: return "double";
        case token::tokenType::WORD: return "word";
        case token::tokenType::DIRECTIVE: return "directive";
        case token::tokenType::STRING: return "string";
        case token::tokenType::EXPRESSION: return "expression";
        case token::tokenType::VARIABLE: return "variable";
        case token::tokenType::VERBATIM: return "verbatim";
        case token::tokenType::CHAR_DATA: return "char_data";
        case token::tokenType::COMPOUND: return "compound";
        case token::tokenType::ERROR: return "error";

        default:
            break;
    }

    return "unknown(" + std::to_string(int(tokType)) + ")";
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, token& tok)
{
    tok.reset();
    return is.read(tok);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token& tok)
{
    switch (tok.type_)
    {
        case token::tokenType::UNDEFINED:
            os << "UNDEFINED";
            WarningInFunction
                << "Undefined token" << endl;
        break;

        case token::tokenType::FLAG:
            // Swallow the flag
        break;

        case token::tokenType::PUNCTUATION:
            os << tok.data_.punctuationVal;
        break;

        case token::tokenType::BOOL:
        case token::tokenType::LABEL:
            os << tok.data_.labelVal;
        break;

        case token::tokenType::FLOAT:
            os << tok.data_.floatVal;
        break;

        case token::tokenType::DOUBLE:
            os << tok.data_.doubleVal;
        break;

        // Possibly different behaviour for serial/parallel streams:
        // preserve types
        case token::tokenType::DIRECTIVE:
        case token::tokenType::EXPRESSION:
        case token::tokenType::VARIABLE:
        case token::tokenType::VERBATIM:
        case token::tokenType::CHAR_DATA:
            os.write(tok);
        break;

        case token::tokenType::WORD:
            os << *tok.data_.wordPtr;
        break;

        case token::tokenType::STRING:
            os << *tok.data_.stringPtr;
        break;

        case token::tokenType::COMPOUND:
            os << *tok.data_.compoundPtr;
        break;

        case token::tokenType::ERROR:
            os << "ERROR";
            WarningInFunction
                << "Error token" << endl;
        break;

        default:
            os << "UNKNOWN";
            SeriousErrorInFunction
                << "Unknown token" << endl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


ostream& Foam::operator<<(ostream& os, const token::punctuationToken& pt)
{
    return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::punctuationToken& pt)
{
    return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::compound& ct)
{
    os  << ct.type();
    if (!ct.pending())
    {
        // Do not write compound content if tagged as 'pending-read'
        os  << token::SPACE;
        ct.write(os);
    }

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ostream& Foam::operator<<
(
    ostream& os,
    const InfoProxy<token>& iproxy
)
{
    return printTokenInfo(os, *iproxy);
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<token>& iproxy
)
{
    return printTokenInfo(os, *iproxy);
}


// ************************************************************************* //
