/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "cellCentreSet.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "volFields.H"
#include "globalIndex.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellCentreSet, 0);
    addToRunTimeSelectionTable(sampledSet, cellCentreSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellCentreSet::genSamples()
{
    const label len = mesh().nCells();

    const globalIndex globalSampleNumbers(len);

    const auto& cellCentres =
        refCast<const fvMesh>(mesh()).C().primitiveField();

    labelList selectedCells = identity(len);
    List<point> selectedPoints;

    if (bounds_.empty())
    {
        selectedPoints = cellCentres;
    }
    else
    {
        label count = 0;
        for (label celli=0; celli < len; ++celli)
        {
            if (bounds_.contains(cellCentres[celli]))
            {
                selectedCells[count++] = celli;
            }
        }

        selectedCells.resize(count);
        selectedPoints = UIndirectList<point>(cellCentres, selectedCells);
    }

    labelList samplingFaces(selectedCells.size(), -1);
    labelList samplingSegments(selectedCells.size(), Zero);
    scalarList samplingCurveDist(selectedCells.size());

    forAll(selectedCells, i)
    {
        samplingCurveDist[i] = globalSampleNumbers.toGlobal(selectedCells[i]);
    }

    // Move into *this
    setSamples
    (
        std::move(selectedPoints),
        std::move(selectedCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellCentreSet::cellCentreSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const boundBox& bbox
)
:
    sampledSet(name, mesh, searchEngine, axis),
    bounds_(bbox)
{
    genSamples();
}


Foam::cellCentreSet::cellCentreSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet
    (
        name,
        mesh,
        searchEngine,
        dict.getOrDefault<word>("axis", "xyz")
    ),
    bounds_(dict.getOrDefault("bounds", boundBox::null()))
{
    genSamples();
}


// ************************************************************************* //
