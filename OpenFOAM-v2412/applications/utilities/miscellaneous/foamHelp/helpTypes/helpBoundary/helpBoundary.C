/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

#include "helpBoundary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace helpTypes
    {
        defineTypeNameAndDebug(helpBoundary, 0);
        addNamedToRunTimeSelectionTable
        (
            helpType,
            helpBoundary,
            dictionary,
            boundary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpTypes::helpBoundary::helpBoundary()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpTypes::helpBoundary::~helpBoundary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpTypes::helpBoundary::init()
{
    helpType::init();

    argList::addOption
    (
        "field",
        "word",
        "List available conditions for field"
    );
    argList::addBoolOption
    (
        "constraint",
        "List constraint patches"
    );
    argList::addBoolOption
    (
        "fixedValue",
        "List fixed value patches (use with -field option)"
    );
}


void Foam::helpTypes::helpBoundary::execute
(
    const argList& args,
    const fvMesh& mesh
)
{
    setEnv("FOAM_ABORT", "", true);

    word condition;
    word fieldName;

    if (args.readIfPresent("browse", condition))
    {
        // TODO: strip scoping info if present?
        // e.g. conditions with leading "compressible::" will not be found
        // ".*[fF]vPatchField.*" + className + ".*"
        displayDoc(condition, ".*[fF]vPatchField.*", false, "H");
    }
    else if (args.found("constraint"))
    {
        wordHashSet constraintTypes(fvPatch::constraintTypes());
        Info<< "Constraint types:" << nl;
        for (const word& cType : constraintTypes)
        {
            Info<< "    " << cType << nl;
        }
        Info<< endl;
    }
    else if (args.readIfPresent("field", fieldName))
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check for any type of volField
        if (fieldHeader.typeHeaderOk<volScalarField>(false))
        {
            if (args.found("fixedValue"))
            {
                fixedValueFieldConditions<scalar>(fieldHeader);
                fixedValueFieldConditions<vector>(fieldHeader);
                fixedValueFieldConditions<sphericalTensor>(fieldHeader);
                fixedValueFieldConditions<symmTensor>(fieldHeader);
                fixedValueFieldConditions<tensor>(fieldHeader);
            }
            else
            {
                (void)fieldConditions<scalar>(fieldHeader, true);
                (void)fieldConditions<vector>(fieldHeader, true);
                (void)fieldConditions<sphericalTensor>(fieldHeader, true);
                (void)fieldConditions<symmTensor>(fieldHeader, true);
                (void)fieldConditions<tensor>(fieldHeader, true);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unable to read field " << fieldName << exit(FatalError);
        }
    }
    else if (args.readIfPresent("fixedValue", fieldName))
    {
        FatalErrorInFunction
            << "-field option must be specified when using the -fixedValue "
            << "option" << exit(FatalError);
    }
    else
    {
        // TODO: strip scoping info if present?
        // e.g. conditions with leading "compressible::" will not be found
        // ".*[fF]vPatchField.*" + className + ".*"
        displayDocOptions(".*[fF]vPatchField.*", false, "H");
    }
}


// ************************************************************************* //
