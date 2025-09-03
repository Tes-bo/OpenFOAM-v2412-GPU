/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "fieldToCell.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "Time.H"
#include "IFstream.H"
#include "fieldDictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, fieldToCell, word);
    addToRunTimeSelectionTable(topoSetSource, fieldToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, fieldToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, fieldToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        fieldToCell,
        word,
        field
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        fieldToCell,
        istream,
        field
    );
}


Foam::topoSetSource::addToUsageTable Foam::fieldToCell::usage_
(
    fieldToCell::typeName,
    "\n    Usage: fieldToCell field min max\n\n"
    "    Select all cells with field value >= min and <= max\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldToCell::applyToSet
(
    const topoSetSource::setAction action,
    const scalarField& field,
    topoSet& set
) const
{
    if (verbose_)
    {
        Info << "    Field min:" << min(field) << " max:" << max(field) << nl;
    }

    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all cells with value of field " << fieldName_
                << " within range " << min_ << ".." << max_ << endl;
        }

        forAll(field, celli)
        {
            if (field[celli] >= min_ && field[celli] <= max_)
            {
                set.set(celli);
            }
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all cells with value of field " << fieldName_
                << " within range " << min_ << ".." << max_ << endl;
        }

        forAll(field, celli)
        {
            if (field[celli] >= min_ && field[celli] <= max_)
            {
                set.unset(celli);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldToCell::fieldToCell
(
    const polyMesh& mesh,
    const word& fieldName,
    const scalar min,
    const scalar max
)
:
    topoSetCellSource(mesh),
    fieldName_(fieldName),
    min_(min),
    max_(max)
{
    if (min_ > max_)
    {
        WarningInFunction
            << "Input min value = " << min_ << " is larger than "
            << "input max value = " << max_ << " for field = " << fieldName_
            << endl;
    }
}


Foam::fieldToCell::fieldToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetCellSource(mesh, dict),
    fieldName_(dict.get<word>("field")),
    min_(dict.get<scalar>("min")),
    max_(dict.get<scalar>("max"))
{}


Foam::fieldToCell::fieldToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetCellSource(mesh),
    fieldName_(checkIs(is)),
    min_(readScalar(checkIs(is))),
    max_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    // Try to load field
    IOobject fieldObject
    (
        fieldName_,
        mesh().time().timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    // Note: should check for volScalarField but that introduces dependency
    //       on volMesh so just use another type with processor-local scope
    if (!fieldObject.typeHeaderOk<labelIOList>(false))
    {
        WarningInFunction
            << "Cannot read field " << fieldName_
            << " from time " << mesh().time().timeName() << endl;
    }
    else if ("volScalarField" == fieldObject.headerClassName())
    {
        // Note: cannot use volScalarField::typeName since that would
        //    introduce linkage problems (finiteVolume needs meshTools)

        IFstream str(fieldObject.typeFilePath<labelIOList>());

        // Read as dictionary
        fieldDictionary fieldDict(fieldObject, fieldObject.headerClassName());

        scalarField internalVals("internalField", fieldDict, mesh().nCells());

        applyToSet(action, internalVals, set);
    }
    else if ("volVectorField" == fieldObject.headerClassName())
    {
        // Note: cannot use volVectorField::typeName since that would
        //    introduce linkage problems (finiteVolume needs meshTools)

        IFstream str(fieldObject.typeFilePath<labelIOList>());

        // Read as dictionary
        fieldDictionary fieldDict(fieldObject, fieldObject.headerClassName());

        vectorField internalVals("internalField", fieldDict, mesh().nCells());

        applyToSet(action, mag(internalVals), set);
    }
    else
    {
        WarningInFunction
            << "Cannot handle fields of type " << fieldObject.headerClassName()
            << endl;
    }
}


// ************************************************************************* //
