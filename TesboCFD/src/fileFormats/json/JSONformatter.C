/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "JSONformatter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(JSONformatter, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::JSONformatter::writeToken(const token& t)
{
    bool ok = true;
    switch (t.type())
    {
        case token::tokenType::BOOL:
            write(t.boolToken());
            break;

        case token::tokenType::LABEL:
            write(t.labelToken());
            break;

        case token::tokenType::FLOAT:
        case token::tokenType::DOUBLE:
            write(t.scalarToken());
            break;

        case token::tokenType::WORD:
        case token::tokenType::DIRECTIVE:
            write(t.wordToken());
            break;

        case token::tokenType::STRING:
        case token::tokenType::EXPRESSION:
        case token::tokenType::VARIABLE:
        case token::tokenType::VERBATIM:
        case token::tokenType::CHAR_DATA:
            write(t.stringToken());
            break;

        default:
            DebugInfo
                << "Problem converting token to JSON:" << nl
                << "    " << Foam::name(int(t.type()))
                << "    - treating as null" << endl;

            os_  << "null";

            ok = false;
            break;
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JSONformatter::JSONformatter(Ostream& os)
:
    os_(os)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::JSONformatter::writeKeyword
(
    const keyType& keyword
)
{
    return os_.writeQuoted(keyword);
}


Foam::Ostream& Foam::JSONformatter::writeDict(const dictionary& dict)
{
    if (dict.empty())
    {
        os_  << "{}";
        return os_;
    }

    os_ << '{' << nl << incrIndent;

    const auto openBrace = [](const token& t){
        return t == token::BEGIN_LIST || t == token::BEGIN_SQR;
    };
    const auto closeBrace = [](const token& t){
        return t == token::END_LIST || t == token::END_SQR;
    };

    label entryi = 0;
    for (const entry& e : dict)
    {
        if (entryi)
        {
            os_  << ',' << nl;
        }

        const word& keyword = e.keyword();

        os_  << indent;
        os_.writeQuoted(keyword) << " : ";

        if (e.isDict())
        {
            writeDict(e.dict());

            ++entryi;

            continue;
        }

        const auto& tokens = e.stream();

        if (tokens.empty()) continue; // error?

        if (tokens.size() == 1)
        {
            writeToken(tokens[0]);
        }
        else
        {
            label offset = 0;
            if (tokens[0].isLabel() && openBrace(tokens[1]))
            {
                offset = 1; // reading 'size (value0 value1 ... valueN)
            }

            const token& t = tokens[offset];
            if (openBrace(t))
            {
                // Assume array-type
                label i = 0;
                for (label tokeni=offset; tokeni<tokens.size(); ++tokeni)
                {
                    const token& tk = tokens[tokeni];

                    if (openBrace(tk))
                    {
                        if (i) os_  << ',';
                        os_  << '[';
                        i = 0;
                    }
                    else if (closeBrace(tk))
                    {
                        os_  << ']';
                    }
                    else
                    {
                        if (i) os_  << ',';
                        if (writeToken(tk)) ++i;

                        if (i % 10 == 0) os_  << nl;
                    }
                }
            }
            else
            {
                // Unknown type - convert to string representation
                os_  << '"';
                forAll(tokens, tokeni)
                {
                    if (tokeni) os_  << ' ';
                    os_  << tokens[tokeni];
                }
                os_  << '"';
            }
        }

        ++entryi;
    }

    os_  << nl << decrIndent << indent << '}';

    return os_;
}


Foam::Ostream& Foam::JSONformatter::write(const bool val)
{
    os_  << (val ? "true" : "false");
    return os_;
}


Foam::Ostream& Foam::JSONformatter::write(const int32_t val)
{
    os_  << val;
    return os_;
}


Foam::Ostream& Foam::JSONformatter::write(const int64_t val)
{
    os_  << val;
    return os_;
}


Foam::Ostream& Foam::JSONformatter::write(const float val)
{
    os_  << val;
    return os_;
}


Foam::Ostream& Foam::JSONformatter::write(const double val)
{
    os_  << val;
    return os_;
}


Foam::Ostream& Foam::JSONformatter::write(const word& str)
{
    return os_.writeQuoted(str);
}


Foam::Ostream& Foam::JSONformatter::write(const std::string& str)
{
    return os_.writeQuoted(str);
}


Foam::Ostream& Foam::JSONformatter::write(const char* c)
{
    os_  << c;
    return os_;
}


// ************************************************************************* //
