/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
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

\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::expressions::exprString::exprString
(
    const std::string& s,
    bool doCheck
)
:
    string(s)
{
    #ifdef FULLDEBUG
    if (doCheck)
    {
        (void)valid();
    }
    #endif
}


inline Foam::expressions::exprString::exprString
(
    std::string&& s,
    bool doCheck
)
:
    string(std::move(s))
{
    #ifdef FULLDEBUG
    if (doCheck)
    {
        (void)valid();
    }
    #endif
}


inline Foam::expressions::exprString::exprString
(
    const char* s,
    bool doCheck
)
:
    string(s)
{
    #ifdef FULLDEBUG
    if (doCheck)
    {
        (void)valid();
    }
    #endif
}


inline Foam::expressions::exprString::exprString
(
    const word& entryName,
    const dictionary& dict,
    const bool mandatory
)
{
    if (mandatory)
    {
        readEntry(entryName, dict);
    }
    else
    {
        readIfPresent(entryName, dict);
    }
}


inline Foam::expressions::exprString::exprString
(
    Istream& is,
    const dictionary& dict,
    const bool removeComments
)
:
    string(is)
{
    expand(dict, removeComments);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::expressions::exprString
Foam::expressions::exprString::toExpr(const std::string& str)
{
    exprString expr;

    expr.string::operator=(str);

    return expr;
}


inline Foam::expressions::exprString
Foam::expressions::exprString::toExpr(std::string&& str)
{
    exprString expr;

    expr.string::operator=(std::move(str));

    return expr;
}


inline Foam::expressions::exprString
Foam::expressions::exprString::toExpr
(
    const std::string& str,
    const dictionary& dict,
    const bool stripComments
)
{
    // Copy - no validate
    exprString expr(str, false);

    expr.expand(dict, stripComments);

    return expr;
}


inline Foam::expressions::exprString
Foam::expressions::exprString::toExpr
(
    std::string&& str,
    const dictionary& dict,
    const bool stripComments
)
{
    // Move - no validate
    exprString expr(std::move(str), false);

    expr.expand(dict, stripComments);

    return expr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::expressions::exprString::valid() const
{
    const bool ok = !contains('$');

    #ifdef FULLDEBUG
    if (!ok)
    {
        FatalErrorInFunction
            << "Unexpanded '$' in " << *this << nl
            << exit(FatalError);
    }
    #endif

    return ok;
}


inline bool Foam::expressions::exprString::readIfPresent
(
    const word& key,
    const dictionary& dict
)
{
    return readEntry(key, dict, IOobjectOption::READ_IF_PRESENT);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::expressions::exprString&
Foam::expressions::exprString::operator=(const char* str)
{
    string::operator=(str);

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


inline Foam::expressions::exprString&
Foam::expressions::exprString::operator=(const std::string& str)
{
    string::operator=(str);

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


inline Foam::expressions::exprString&
Foam::expressions::exprString::operator=(std::string&& str)
{
    string::operator=(std::move(str));

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


// ************************************************************************* //
