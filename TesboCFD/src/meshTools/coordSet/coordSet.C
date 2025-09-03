/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "coordSet.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::coordSet::coordFormat
>
Foam::coordSet::coordFormatNames
({
    { coordFormat::X, "x" },
    { coordFormat::Y, "y" },
    { coordFormat::Z, "z" },
    { coordFormat::RADIUS, "radius" },
    { coordFormat::DISTANCE, "distance" },
    { coordFormat::XYZ, "xyz" },
/// { coordFormat::DEFAULT, "default" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordSet::checkDimensions() const
{
    if (points().size() != distance().size())
    {
        FatalErrorInFunction
            << "Size not equal :" << nl
            << "    points:" << points().size()
            << "    distance:" << distance().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Foam::coordSet::coordSet()
// :
//     pointField(),
//     name_(),
//     distance_(),
//     axis_(coordFormat::DEFAULT)
// {}


Foam::coordSet::coordSet
(
    const word& name,
    const coordFormat axisType
)
:
    pointField(),
    name_(name),
    distance_(),
    axis_(axisType)
{}


Foam::coordSet::coordSet
(
    const word& name,
    const word& axis
)
:
    coordSet(name, coordFormatNames.get(axis))
{}


Foam::coordSet::coordSet
(
    const word& name,
    const word& axis,
    const List<point>& points,
    const scalarList& dist
)
:
    pointField(points),
    name_(name),
    distance_(dist),
    axis_(coordFormatNames[axis])
{
    checkDimensions();
}


Foam::coordSet::coordSet
(
    const word& name,
    const word& axis,
    List<point>&& points,
    scalarList&& dist
)
:
    pointField(std::move(points)),
    name_(name),
    distance_(std::move(dist)),
    axis_(coordFormatNames.get(axis))
{
    checkDimensions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSet::hasVectorAxis() const noexcept
{
    return axis_ == coordFormat::XYZ;
}


Foam::scalar Foam::coordSet::scalarCoord(const label index) const
{
    switch (axis_)
    {
        case coordFormat::X:
        {
            return points()[index].x();
        }
        case coordFormat::Y:
        {
            return points()[index].y();
        }
        case coordFormat::Z:
        {
            return points()[index].z();
        }
        case coordFormat::RADIUS:
        {
            return mag(points()[index]);
        }
        case coordFormat::DISTANCE:
        {
            // Note: the distance will unset it constructed from
            // 'name' and 'axis' only

            if (distance().empty())
            {
                FatalErrorInFunction
                    << "Axis type '" << coordFormatNames[axis_]
                    << "' requested but curve distance has not been set"
                    << abort(FatalError);
            }
            return distance()[index];
        }
        default:
        {
            FatalErrorInFunction
                << "Illegal axis specification '" << coordFormatNames[axis_]
                << "' for sampling " << name_
                << exit(FatalError);
        }
    }

    return 0;
}


const Foam::point& Foam::coordSet::vectorCoord(const label index) const
{
    return points()[index];
}


Foam::Ostream& Foam::coordSet::write(Ostream& os) const
{
    os  << "name:" << name_ << " axis:" << coordFormatNames[axis_] << nl
        << nl
        << "\t(coord)" << nl;

    for (const point& p : *this)
    {
        os  << '\t' << p << nl;
    }

    return os;
}


Foam::autoPtr<Foam::coordSet> Foam::coordSet::gatherSort
(
    labelList& sortOrder
) const
{
    // Combine sampleSet from processors. Sort by distance.
    // Return ordering in indexSet.
    // Note: only master results are valid

    List<point> allPoints(globalIndex::gatherOp(points()));
    List<scalar> allDistance(globalIndex::gatherOp(distance()));

    if (Pstream::master() && allDistance.empty())
    {
        WarningInFunction
            << "Gathered empty coordSet: " << name() << endl;
    }

    // Sort according to distance
    Foam::sortedOrder(allDistance, sortOrder);  // uses stable sort

    // Repopulate gathered points/distances in the correct order
    allPoints = List<point>(allPoints, sortOrder);
    allDistance = List<scalar>(allDistance, sortOrder);

    return autoPtr<coordSet>::New
    (
        name(),
        axis(),
        std::move(allPoints),
        std::move(allDistance)
    );
}


// ************************************************************************* //
