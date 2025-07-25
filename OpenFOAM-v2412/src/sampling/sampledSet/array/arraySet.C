/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018,2022 OpenCFD Ltd.
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

#include "arraySet.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(arraySet, 0);
    addToRunTimeSelectionTable(sampledSet, arraySet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::arraySet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    const meshSearch& queryMesh = searchEngine();

    if (cmptMin(pointsDensity_) < 1)
    {
        FatalErrorInFunction << "Illegal specification for pointsDensity "
            << pointsDensity_
            << ". It specifies the number of points per coordinate axis"
            << " so should contain positive integers only."
            << exit(FatalError);
    }

    const scalar dx = spanBox_.x()/(pointsDensity_.x() + 1);
    const scalar dy = spanBox_.y()/(pointsDensity_.y() + 1);
    const scalar dz = spanBox_.z()/(pointsDensity_.z() + 1);

    label sampleI(0);

    for (label k=1; k<=pointsDensity_.z(); ++k)
    {
        for (label j=1; j<=pointsDensity_.y(); ++j)
        {
            for (label i=1; i<=pointsDensity_.x(); ++i)
            {
                // Local Cartesian
                point pt(i*dx, j*dy, k*dz);

                // Global Cartesian
                pt = csys_.globalPosition(pt);

                const label celli = queryMesh.findCell(pt);

                if (celli != -1)
                {
                    samplingPts.append(pt);
                    samplingCells.append(celli);
                    samplingFaces.append(-1);
                    samplingSegments.append(0);
                    samplingCurveDist.append(1.0 * sampleI);
                }
            }
        }
    }
}


void Foam::arraySet::genSamples()
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

Foam::arraySet::arraySet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const coordSystem::cartesian& csys,
    const Vector<label>& pointsDensity,
    const Vector<scalar>& spanBox
)
:
    sampledSet(name, mesh, searchEngine, axis),
    csys_(csys),
    pointsDensity_(pointsDensity),
    spanBox_(spanBox)
{
    genSamples();
}


Foam::arraySet::arraySet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    csys_(dict),  // Note: no indirect cs with this constructor
    pointsDensity_(dict.get<labelVector>("pointsDensity")),
    spanBox_(dict.get<vector>("spanBox"))
{
    genSamples();
}


// ************************************************************************* //
