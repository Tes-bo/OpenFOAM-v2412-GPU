/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "circleSet.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(circleSet, 0);
    addToRunTimeSelectionTable(sampledSet, circleSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::circleSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    // Set start point
    label celli = searchEngine().findCell(startPoint_);
    if (celli != -1)
    {
        samplingPts.append(startPoint_);
        samplingCells.append(celli);
        samplingFaces.append(-1);
        samplingSegments.append(0);
        samplingCurveDist.append(0.0);
    }
    else
    {
        WarningInFunction
            << "Unable to find cell at point id " << 0
            << " at location " << startPoint_ << endl;
    }

    // Add remaining points
    const scalar alpha = degToRad(dTheta_);
    const scalar sinAlpha = sin(alpha);
    const scalar cosAlpha = cos(alpha);

    // First axis
    vector axis1 = startPoint_ - origin_;
    const scalar radius = mag(axis1);

    if (mag(axis1 & circleAxis_) > SMALL)
    {
        WarningInFunction
            << "Vector defined by (startPoint - origin) not orthogonal to "
            << "circleAxis:" << nl
            << "    startPoint - origin = " << axis1 << nl
            << "    circleAxis          = " << circleAxis_ << nl
            << endl;
    }

    axis1 /= mag(axis1);

    scalar theta = dTheta_;
    label nPoint = 1;
    while (theta < 360)
    {
        axis1 = normalised(axis1*cosAlpha + (axis1^circleAxis_)*sinAlpha);
        point pt = origin_ + radius*axis1;

        label celli = searchEngine().findCell(pt);

        if (celli != -1)
        {
            samplingPts.append(pt);
            samplingCells.append(celli);
            samplingFaces.append(-1);
            samplingSegments.append(nPoint);
            samplingCurveDist.append(radius*degToRad(theta));

            ++nPoint;
        }
        else
        {
            WarningInFunction
                << "Unable to find cell at point id " << nPoint
                << " at location " << pt << endl;
        }
        theta += dTheta_;
    }
}


void Foam::circleSet::genSamples()
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

Foam::circleSet::circleSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const point& origin,
    const vector& circleAxis,
    const point& startPoint,
    const scalar dTheta
)
:
    sampledSet(name, mesh, searchEngine, axis),
    origin_(origin),
    circleAxis_(circleAxis),
    startPoint_(startPoint),
    dTheta_(dTheta)
{
    genSamples();
}


Foam::circleSet::circleSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    origin_(dict.get<point>("origin")),
    circleAxis_(normalised(dict.get<vector>("circleAxis"))),
    startPoint_(dict.get<point>("startPoint")),
    dTheta_(dict.get<scalar>("dTheta"))
{
    genSamples();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::circleSet::getRefPoint(const List<point>& pts) const
{
    return startPoint_;
}


// ************************************************************************* //
