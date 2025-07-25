/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "ddt2.H"
#include "stringOps.H"
#include "volFields.H"
#include "dictionary.H"
#include "wordRes.H"
#include "steadyStateDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ddt2, 0);
    addToRunTimeSelectionTable(functionObject, ddt2, dictionary);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Check that string contains the appropriate substitution token(s)
static bool checkFormatName(const word& str)
{
    if (!str.contains("@@"))
    {
        WarningInFunction
            << "Bad result naming (no '@@' token found)."
            << nl << endl;

        return false;
    }
    else if (str == "@@")
    {
        WarningInFunction
            << "Bad result naming (only a '@@' token found)."
            << nl << endl;

        return false;
    }

    return true;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ddt2::accept(const word& fieldName) const
{
    // check input vs possible result names
    // to avoid circular calculations
    return !denyField_.match(fieldName);
}


int Foam::functionObjects::ddt2::process(const word& fieldName)
{
    if (!accept(fieldName))
    {
        return -1;
    }

    int state = 0;

    apply<volScalarField>(fieldName, state);
    apply<volVectorField>(fieldName, state);

    return state;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ddt2::ddt2
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    selectFields_(),
    resultName_(),
    denyField_(),
    results_(),
    mag_(dict.getOrDefault("mag", false))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ddt2::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (word(mesh_.ddtScheme("default")) == "steadyState")
    {
        WarningInFunction
            << typeName << " function object not appropriate for steady-state"
            << endl;
        return false;
    }

    dict.readEntry("fields", selectFields_);
    selectFields_.uniq();

    Info<< type() << " fields: " << selectFields_ << nl;

    resultName_ = dict.getOrDefault<word>
    (
        "result",
        ( mag_ ? "mag(ddt(@@))" : "magSqr(ddt(@@))" )
    );

    // Expect '@@' token for result, unless a single (non-regex) source field
    if
    (
        (selectFields_.size() == 1 && selectFields_.first().isLiteral())
     || checkFormatName(resultName_)
    )
    {
        denyField_.set
        (
            stringOps::quotemeta(resultName_, regExp::meta())
            .replace("@@", "(.+)")
        );

        return true;
    }

    denyField_.clear();
    return false;
}


bool Foam::functionObjects::ddt2::execute()
{
    results_.clear();

    wordHashSet candidates(mesh_.names<regIOobject>(selectFields_));

    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // Check exact matches first
    for (const wordRe& select : selectFields_)
    {
        if (select.isLiteral())
        {
            const word& fieldName = select;

            if (!candidates.erase(fieldName))
            {
                missing.append(fieldName);
            }
            else if (process(fieldName) < 1)
            {
                ignored.append(fieldName);
            }
        }
    }

    for (const word& fieldName : candidates)
    {
        process(fieldName);
    }

    if (missing.size())
    {
        WarningInFunction
            << "Missing field " << missing << endl;
    }
    if (ignored.size())
    {
        WarningInFunction
            << "Unprocessed field " << ignored << endl;
    }

    return true;
}


bool Foam::functionObjects::ddt2::write()
{
    if (results_.size())
    {
        Log << type() << ' ' << name() << " write:" << endl;
    }

    // Consistent output order
    const wordList outputList = results_.sortedToc();
    for (const word& fieldName : outputList)
    {
        if (foundObject<regIOobject>(fieldName))
        {
            const regIOobject& io = lookupObject<regIOobject>(fieldName);

            Log << "    " << fieldName << endl;

            io.write();
        }
    }

    return true;
}


// ************************************************************************* //
