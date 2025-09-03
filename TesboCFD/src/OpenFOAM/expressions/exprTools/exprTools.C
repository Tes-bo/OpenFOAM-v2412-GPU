/*---------------------------------------------------------------------------*\
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

Note
    Ideas based on swak4Foam driver code (2010-2018)
    from Bernhard Gschaider <bgschaid@hfd-research.com>

\*---------------------------------------------------------------------------*/

#include "exprTools.H"
#include "stringOps.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//! \cond file_scope

// Maximum depth for recursive variable names
static constexpr label maxRecursionDepth_ = 100;


static List<expressions::exprString> expandExprStrings
(
    const UList<string>& inputs,
    const dictionary& dict,
    bool mandatory,
    label recursionDepth
)
{
    ///Info<< "::expandExprStrings " << inputs << endl;

    DynamicList<expressions::exprString> result;

    for (const string& input : inputs)
    {
        // Allow inline list of semicolon-separated variables
        const auto varExpressions = stringOps::split<string>(input, ';');

        for (const auto& subMatch : varExpressions)
        {
            string varExpr(stringOps::trim(subMatch.str()));

            if (varExpr.empty())
            {
                continue;
            }

            ///Info<< "Checking " << varExpr << endl;

            // Expand #otherVariable as dictionary lookup
            if (varExpr[0] == '#')
            {
                ///Info<< "Expand: " << varExpr << endl;

                List<expressions::exprString> expansions
                (
                    exprTools::getList
                    (
                        varExpr.substr(1),
                        dict,
                        mandatory,
                        recursionDepth
                    )
                );

                ///Info<< "Expanded " << varExpr.c_string() << ": "
                ///    << expansions << nl;

                result.reserve(result.size() + expansions.size());
                for (expressions::exprString& str : expansions)
                {
                    result.append(std::move(str));
                }
            }
            else
            {
                result.append
                (
                    expressions::exprString::toExpr(std::move(varExpr), dict)
                );
            }
        }
    }

    result.shrink();
    return result;
}

//! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::List<Foam::expressions::exprString>
Foam::exprTools::getList
(
    const word& keyword,
    const dictionary& dict,
    bool mandatory,
    label recursionDepth
)
{
    List<expressions::exprString> result;

    // Catch empty keyword as a no-op (eg, when called recursively)
    if (keyword.empty())
    {
        return result;
    }

    const entry* eptr = dict.findEntry(keyword, keyType::LITERAL_RECURSIVE);

    if (!eptr)
    {
        if (mandatory)
        {
            FatalIOErrorInFunction(dict)
                << "Missing mandatory entry: " << keyword << nl << nl
                << exit(FatalIOError);
        }

        return result;
    }

    if (++recursionDepth > maxRecursionDepth_)
    {
        FatalIOErrorInFunction(dict)
            << "Exceeded recursion depth (" << maxRecursionDepth_
            << ") while reading list " << keyword << nl
            << "Likely caused by circular referencing" << nl
            << exit(FatalIOError);
    }


    ITstream& is = eptr->stream();
    token tok(is);

    List<string> list;

    if (tok.isLabel() || tok.isPunctuation(token::BEGIN_LIST))
    {
        // A list of strings
        is.rewind();
        is >> list;
    }
    else if (tok.isString())
    {
        // A single string
        list.resize(1);
        list[0] = tok.stringToken();
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << " Entry '"<< keyword
            << "' not a string or list of strings" << nl
            << exit(FatalIOError);

        return result;
    }

    // Check for excess tokens
    dict.checkITstream(is, keyword);

    // Expand List<string> to List<expressions::exprString>
    return expandExprStrings(list, dict, mandatory, recursionDepth);
}


// ************************************************************************* //
