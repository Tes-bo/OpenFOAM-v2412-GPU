/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "abaqusMeshSet.H"
#include "stringOps.H"
#include "Fstream.H"
#include "SpanStream.H"
#include "meshSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(abaqusMeshSet, 0);
    addToRunTimeSelectionTable(sampledSet, abaqusMeshSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::abaqusMeshSet::readCoord(ISstream& is, vector& coord) const
{
    string buffer;

    do
    {
        //buffer.clear();

        is.getLine(buffer);

        const auto elems = buffer.find("*ELEMENT");

        if (elems != std::string::npos)
        {
            buffer.clear();
            break;
        }

        // Trim out any '*' comments
        const auto pos = buffer.find('*');
        if (pos != std::string::npos)
        {
            buffer.erase(pos);
        }
        stringOps::inplaceTrimRight(buffer);
    }
    while (buffer.empty() && is.good());

    if (buffer.empty())
    {
        return false;
    }

    const auto strings = stringOps::split(buffer, ',');

    if (strings.size() != 4)
    {
        FatalErrorInFunction
            << "Read error: expected format int, float, float, float"
            << " but read buffer " << buffer
            << exit(FatalError);
    }

    for (int i = 0; i <= 2; ++i)
    {
        // Swallow i=0 for node label
        const auto& s = strings[i+1].str();
        ISpanStream buf(s.data(), s.length());
        buf >> coord[i];
    }

    coord *= scale_;

    return true;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::abaqusMeshSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    DebugInfo
        << "abaqusMeshSet : sampling " << sampleCoords_.size() << " points"
        << endl;

    List<bool> found(sampleCoords_.size(), false);

    forAll(sampleCoords_, samplei)
    {
        const vector& pt = sampleCoords_[samplei];

        label celli = searchEngine().findCell(pt);
        if (celli != -1)
        {
            found[samplei] = true;

            samplingPts.append(pt);
            samplingCells.append(celli);
            samplingFaces.append(-1);
            samplingSegments.append(0);
            samplingCurveDist.append(1.0*samplei);
        }
    }

    Pstream::listCombineAllGather(found, orEqOp<bool>());

    DynamicList<label> lost;
    forAll(found, samplei)
    {
        if (!found[samplei]) lost.append(samplei);
    }

    const label nFound = sampleCoords_.size() - lost.size();
    label nOutOfBounds = 0;

    if (Pstream::parRun())
    {
        typedef Tuple2<label, Tuple2<label, scalar>> procCellDist;
        List<procCellDist> pcd(lost.size(), procCellDist(-1, {-1, GREAT}));
        forAll(pcd, i)
        {
            const label samplei = lost[i];
            const vector& pt0 = sampleCoords_[samplei];
            const label celli = searchEngine().findNearestCell(pt0);
            const vector& pt1 = mesh().cellCentres()[celli];
            const scalar distSqr = magSqr(pt1 - pt0);
            pcd[i] = procCellDist(Pstream::myProcNo(), {celli, distSqr});
        }

        Pstream::listCombineAllGather
        (
            pcd,
            [](procCellDist& x, const procCellDist& y)
            {
                // Select item with smallest distance
                if (y.second().second() < x.second().second())
                {
                    x = y;
                }
            }
        );

        forAll(pcd, i)
        {
            const label samplei = lost[i];

            if (pcd[i].second().second() < maxDistSqr_)
            {
                if (pcd[i].first() == Pstream::myProcNo())
                {
                    const label celli = pcd[i].second().first();
                    const vector& pt1 = mesh().cellCentres()[celli];

                    samplingPts.append(pt1);
                    samplingCells.append(celli);
                    samplingFaces.append(-1);
                    samplingSegments.append(0);
                    samplingCurveDist.append(1.0*samplei);
                }
            }
            else
            {
                // Insert points that have not been found as null points
                if (Pstream::master())
                {
                    samplingPts.append(sampleCoords_[samplei]);
                    samplingCells.append(-1);
                    samplingFaces.append(-1);
                    samplingSegments.append(0);
                    samplingCurveDist.append(1.0*samplei);
                    ++nOutOfBounds;
                }
            }
        }
    }
    else
    {
        // Serial running

        forAll(lost, i)
        {
            const label samplei = lost[i];
            const vector& pt0 = sampleCoords_[samplei];
            const label celli = searchEngine().findNearestCell(pt0);
            const vector& pt1 = mesh().cellCentres()[celli];
            const scalar distSqr = magSqr(pt1 - pt0);

            if (distSqr < maxDistSqr_)
            {
                samplingPts.append(pt1);
                samplingCells.append(celli);
                samplingFaces.append(-1);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0*samplei);
            }
            else
            {
                // Insert points that have not been found as null points
                samplingPts.append(sampleCoords_[samplei]);
                samplingCells.append(-1);
                samplingFaces.append(-1);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0*samplei);
                ++nOutOfBounds;
            }
        }
    }

    DebugInFunction
        << "Sample size       : " << sampleCoords_.size() << nl
        << "Lost samples      : " << nOutOfBounds << nl
        << "Recovered samples : "
        << (sampleCoords_.size() - nOutOfBounds - nFound) << nl
        << endl;


    if (nOutOfBounds)
    {
        WarningInFunction
            << "Identified " << nOutOfBounds << " out-of-bounds points"
            << endl;
    }
}


void Foam::abaqusMeshSet::genSamples()
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

    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug > 1)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::abaqusMeshSet::abaqusMeshSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    scale_(dict.getOrDefault<scalar>("scale", 1)),
    sampleCoords_(),
    maxDistSqr_(sqr(dict.getOrDefault<scalar>("maxDist", 0)))
{
    if (Pstream::master())
    {
        const fileName inputFile(dict.get<fileName>("file").expand());
        IFstream pointsFile(inputFile);

        if (!pointsFile.good())
        {
            FatalIOErrorInFunction(dict)
                << "Unable to find file " << pointsFile.name()
                << abort(FatalIOError);
        }

        // Read the points file
        DynamicList<point> coords;
        vector c;
        while (readCoord(pointsFile, c))
        {
            coords.append(c);
        }

        sampleCoords_.transfer(coords);
    }

    Pstream::broadcast(sampleCoords_);

    DebugInfo
        << "Number of sample points: " << sampleCoords_.size() << nl
        << "Sample points bounds: " << boundBox(sampleCoords_) << endl;

    genSamples();
}


// ************************************************************************* //
