/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "cloudSet.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, cloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    const meshSearch& queryMesh = searchEngine();

    labelList foundProc(sampleCoords_.size(), -1);
    forAll(sampleCoords_, sampleI)
    {
        label celli = queryMesh.findCell(sampleCoords_[sampleI]);

        if (celli != -1)
        {
            samplingPts.append(sampleCoords_[sampleI]);
            samplingCells.append(celli);
            samplingFaces.append(-1);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0 * sampleI);

            foundProc[sampleI] = Pstream::myProcNo();
        }
    }

    // Check that all have been found
    labelList maxFoundProc(foundProc);
    labelList minFoundProc(foundProc.size(), labelMax);
    forAll(foundProc, i)
    {
        if (foundProc[i] != -1)
        {
            minFoundProc[i] = foundProc[i];
        }
    }
    Pstream::listCombineReduce(minFoundProc, minEqOp<label>());
    Pstream::listCombineReduce(maxFoundProc, maxEqOp<label>());


    DynamicField<point> missingPoints(sampleCoords_.size());

    forAll(sampleCoords_, sampleI)
    {
        if (maxFoundProc[sampleI] == -1)
        {
            // No processor has found the location.
            missingPoints.append(sampleCoords_[sampleI]);
        }
        else if (minFoundProc[sampleI] != maxFoundProc[sampleI])
        {
            WarningInFunction
                << "For sample set " << name()
                << " location " << sampleCoords_[sampleI]
                << " seems to be on multiple domains: "
                << minFoundProc[sampleI] << " and " << maxFoundProc[sampleI]
                << nl
                << "This might happen if the location is on"
                << " a processor patch. Change the location slightly"
                << " to prevent this." << endl;
        }
    }


    if (missingPoints.size() > 0)
    {
        if (missingPoints.size() < 100 || debug)
        {
            WarningInFunction
                << "For sample set " << name()
                << " did not found " << missingPoints.size()
                << " points out of " << sampleCoords_.size()
                << nl
                << "Missing points:" << missingPoints << endl;
        }
        else
        {
            WarningInFunction
                << "For sample set " << name()
                << " did not found " << missingPoints.size()
                << " points out of " << sampleCoords_.size()
                << nl
                << "Print missing points by setting the debug flag"
                << " for " << cloudSet::typeName << endl;
        }
    }
}


void Foam::cloudSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
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

Foam::cloudSet::cloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords)
{
    genSamples();
}


Foam::cloudSet::cloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.get<pointField>("points"))
{
    genSamples();
}


// ************************************************************************* //
