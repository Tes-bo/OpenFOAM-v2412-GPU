/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 Bernhard Gschaider
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

\*---------------------------------------------------------------------------*/

#include "expressionEntry.H"
#include "stringOps.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

defineTypeName(expressionEntry);
defineRunTimeSelectionTable(expressionEntry, empty);

// Various types can be used directly without any changes

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    direct
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    label
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    scalar
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    word
);

} // End namespace exprTools
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
inline static const entry* getVariableOrDie
(
    const word& name,
    const dictionary& dict
)
{
    const entry* eptr = dict.findScoped(name, keyType::LITERAL_RECURSIVE);

    if (!eptr)
    {
        FatalIOErrorInFunction(dict)
            << "No dictionary entry " << name << nl
            << exit(FatalIOError);
    }

    if (eptr->isDict())
    {
        FatalIOErrorInFunction(dict)
            << "Found dictionary " << name << " instead of entry" << nl
            << exit(FatalIOError);
    }

    return eptr;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::exprTools::expressionEntry>
Foam::exprTools::expressionEntry::New
(
    const word& name
)
{
    auto* ctorPtr = emptyConstructorTable(name);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "expressionEntry",
            name,
            *emptyConstructorTablePtr_
        )  << exit(FatalError);
    }

    return autoPtr<expressionEntry>(ctorPtr());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::exprTools::expressionEntry::inplaceExpand
(
    std::string& s,
    const dictionary& dict
)
{
    // This is much like stringOps::inplaceExpand
    constexpr const char sigil = '$';

    // Step 1:
    // Handle $[] special expansions first

    std::string::size_type varBeg = 0;

    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg && s[varBeg-1] == '\\')
        {
            // Escaped character - pass through
            ++varBeg;
            continue;
        }

        if (s[varBeg+1] == '[')
        {
            // An expression pattern with $[...]

            std::string::size_type varEnd = s.find(']', varBeg);
            std::string::size_type delim = 1;

            if (varEnd == std::string::npos)
            {
                // Parsed '$[...' without closing ']' - error
                FatalErrorInFunction
                    << "No correct terminating ']' found in " << s << nl
                    << exit(FatalError);
                break;
            }

            // Look for embedded (type) cast
            word castTo, varName;

            const auto lparen = varBeg+2;
            if (lparen < s.size() && s[lparen] == '(')
            {
                const auto rparen = s.find(')', lparen);

                if (rparen > varEnd)
                {
                    // Handles both "$[( ...]" and "$[( ...])" cases

                    auto& err = FatalErrorInFunction;

                    if (rparen == std::string::npos)
                    {
                        err << "No closing ')' found in ";
                    }
                    else
                    {
                        err << "Closing ')' found outside of";
                    }

                    err << " substring "
                        << s.substr(varBeg, varEnd-varBeg) << nl
                        << exit(FatalError);
                }

                castTo.assign(s.substr(lparen+1, rparen - lparen - 1));
                varName.assign(s.substr(rparen+1, varEnd - rparen - 1));
            }
            else
            {
                varName.assign
                (
                    s.substr(varBeg + 1 + delim, varEnd - varBeg - 2*delim)
                );
            }

            // Likely no spaces there, but for extra safety...
            stringOps::inplaceTrim(varName);

            // Allow recursive plain expansion for the *variable* name.
            // This means "$[(vector) var${index} ]" should work

            // Expand with env=true, empty=true, subDict=false
            stringOps::inplaceExpand(varName, dict, true, true, false);

            // Length of original text to replace (incl. decorators)
            const auto replaceLen = (varEnd - varBeg + 1);

            // Get primitiveEntry with env=false, subDict=false
            const entry* eptr = getVariableOrDie(varName, dict);

            std::string varValue;

            if (castTo.empty())
            {
                // Serialized with spaces - fails for non-primitiveEntry
                ITstream& its = eptr->stream();
                varValue = its.toString();
            }
            else
            {
                varValue = expressionEntry::New(castTo)->toExpr(*eptr);
            }

            s.std::string::replace(varBeg, replaceLen, varValue);
            varBeg += varValue.size();
        }
        else
        {
            ++varBeg;
        }
    }


    // Step 2:
    // Handle all ${}, $var and ${{ ... }} expansions.
    // - this is done second such that $[(vector) xyz] entries will have
    //   been properly expanded by this stage

    // Expand with env=true, empty=true, subDict=false
    stringOps::inplaceExpand(s, dict, true, true, false);
}


Foam::expressions::exprString
Foam::exprTools::expressionEntry::expand
(
    const std::string& orig,
    const dictionary& dict
)
{
    // Copy without validation (use assign)
    expressions::exprString s;
    s.assign(orig);

    inplaceExpand(s, dict);

    return s;
}


// ************************************************************************* //
