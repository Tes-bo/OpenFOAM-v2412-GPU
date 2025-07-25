/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "solverFieldSelection.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solverFieldSelection::solverFieldSelection
(
    const objectRegistry& obr,
    const bool includeComponents
)
:
    fieldSelection(obr, includeComponents)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "Registry must be of type " << fvMesh::typeName
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::solverFieldSelection::updateSelection()
{
    List<fieldInfo> oldSet(std::move(selection_));

    DynamicList<fieldInfo> newSelection(oldSet.size());

    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const dictionary& solverDict = mesh.data().solverPerformanceDict();

    const wordList solvedFieldNames(solverDict.sortedToc());

    for (const fieldInfo& fi : *this)
    {
        for (const word& solvedField : solvedFieldNames)
        {
            if (fi.name().match(solvedField))
            {
                fi.found(true);
                newSelection.emplace_back(wordRe(solvedField), fi.component());
            }
        }
    }

    selection_.transfer(newSelection);

    if (!fieldSelection::checkSelection())
    {
        WarningInFunction
            << "Valid solver fields are: " << solvedFieldNames;
    }

    return selection_ != oldSet;
}


// ************************************************************************* //
