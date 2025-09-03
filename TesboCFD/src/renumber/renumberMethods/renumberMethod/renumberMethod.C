/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "renumberMethod.H"
#include "dlLibraryTable.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(renumberMethod);
    defineRunTimeSelectionTable(renumberMethod, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::wordList Foam::renumberMethod::supportedMethods()
{
    if (dictionaryConstructorTablePtr_)
    {
        return dictionaryConstructorTablePtr_->sortedToc();
    }
    return wordList();
}


Foam::autoPtr<Foam::renumberMethod> Foam::renumberMethod::New
(
    const dictionary& dict
)
{
    const word methodType(dict.get<word>("method"));

    //Info<< "Selecting renumberMethod " << methodType << endl;

    // Load any additional libs (verbose = false)
    dlLibraryTable::libs().open("libs", dict, false);

    auto* ctorPtr = dictionaryConstructorTable(methodType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "renumberMethod",
            methodType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<renumberMethod>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::renumberMethod::renumber(const label nCells) const
{
    NotImplemented;
    return labelList();
}


Foam::labelList Foam::renumberMethod::renumber(const pointField& cc) const
{
    NotImplemented;
    return labelList();
}


Foam::labelList Foam::renumberMethod::renumber
(
    const polyMesh& mesh
) const
{
    // Local mesh connectivity
    CompactListList<label> cellCells;
    globalMeshData::calcCellCells(mesh, cellCells);

    return renumber(cellCells);
}


Foam::labelList Foam::renumberMethod::renumber
(
    const polyMesh& mesh,
    const labelUList& fineToCoarse,
    const pointField& coarsePoints
) const
{
    CompactListList<label> coarseCellCells;
    globalMeshData::calcCellCells
    (
        mesh,
        fineToCoarse,
        coarsePoints.size(),   // nLocalCoarse
        false,  // local only (parallel = false)
        coarseCellCells
    );

    // Renumber based on agglomerated points
    labelList coarseDistribution = renumber(coarseCellCells);

    // From coarse back to fine for original mesh
    return labelList(coarseDistribution, fineToCoarse);
}


// ************************************************************************* //
