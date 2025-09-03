/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "searchableRotatedBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableRotatedBox, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableRotatedBox,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableRotatedBox,
        dict,
        rotatedBox
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableRotatedBox::searchableRotatedBox
(
    const IOobject& io,
    const vector& span,
    const coordSystem::cartesian& csys
)
:
    searchableSurface(io),
    box_
    (
        IOobject
        (
            io.name() + "_box",
            io.instance(),
            io.local(),
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false      // never register
        ),
        treeBoundBox(Zero, span)
    ),
    transform_(csys.origin(), csys.e3(), csys.e1())
{
    points_ = transform_.globalPosition(box_.points());
}


Foam::searchableRotatedBox::searchableRotatedBox
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableRotatedBox
    (
        io,
        dict.get<vector>("span"),
        coordSystem::cartesian
        (
            dict.get<point>("origin"),
            dict.get<vector>("e3"),
            dict.get<vector>("e1")
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableRotatedBox::regions() const
{
    return box_.regions();
}


Foam::tmp<Foam::pointField> Foam::searchableRotatedBox::coordinates() const
{
    return transform_.globalPosition(box_.coordinates());
}


void Foam::searchableRotatedBox::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    box_.boundingSpheres(centres, radiusSqr);
    centres = transform_.globalPosition(centres);
}


Foam::tmp<Foam::pointField> Foam::searchableRotatedBox::points() const
{
    return points_;
}


bool Foam::searchableRotatedBox::overlaps(const boundBox& bb) const
{
    // (from treeDataPrimitivePatch.C)

    // 1. bounding box
    if (!bb.overlaps(bounds()))
    {
        return false;
    }

    // 2. Check if any corner points inside
    if (bb.containsAny(points_))
    {
        return true;
    }

    // 3. Difficult case: all points are outside but connecting edges might
    // go through cube.

    const treeBoundBox treeBb(bb);

    // 3a. my edges through bb faces
    for (const edge& e : treeBoundBox::edges)
    {
        point inter;
        if (treeBb.intersects(points_[e[0]], points_[e[1]], inter))
        {
            return true;
        }
    }

    // 3b. bb edges through my faces

    const pointField bbPoints(bb.points());

    for (const face& f : treeBoundBox::faces)
    {
        const point fc = f.centre(points_);

        for (const edge& e : treeBoundBox::edges)
        {
            pointHit inter = f.intersection
            (
                bbPoints[e[0]],
                bbPoints[e[1]],
                fc,
                points_,
                intersection::HALF_RAY
            );

            if (inter.hit() && inter.distance() <= 1)
            {
                return true;
            }
        }
    }

    return false;
}


Foam::pointIndexHit Foam::searchableRotatedBox::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit boxNearest
    (
        box_.findNearest
        (
            transform_.localPosition(sample),
            nearestDistSqr
        )
    );

    boxNearest.point() = transform_.globalPosition(boxNearest.point());

    return boxNearest;
}


Foam::pointIndexHit Foam::searchableRotatedBox::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    NotImplemented;
    return pointIndexHit();
}


Foam::pointIndexHit Foam::searchableRotatedBox::findLine
(
    const point& start,
    const point& end
) const
{
    pointIndexHit boxHit
    (
        box_.findLine
        (
            transform_.localPosition(start),
            transform_.localPosition(end)
        )
    );

    boxHit.point() = transform_.globalPosition(boxHit.point());

    return boxHit;
}


Foam::pointIndexHit Foam::searchableRotatedBox::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine(start, end);
}


void Foam::searchableRotatedBox::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableRotatedBox::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        info[i] = findLine(start[i], end[i]);
    }
}


void Foam::searchableRotatedBox::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        info[i] = findLineAny(start[i], end[i]);
    }
}


void Foam::searchableRotatedBox::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    // Work array
    DynamicList<pointIndexHit> hits;

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress
    const vectorField dirVec(end-start);
    const scalarField magSqrDirVec(magSqr(dirVec));
    const vectorField smallVec
    (
        Foam::sqrt(SMALL)*dirVec
      + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
    );

    forAll(start, pointI)
    {
        // See if any intersection between pt and end
        pointIndexHit inter = findLine(start[pointI], end[pointI]);

        if (inter.hit())
        {
            hits.clear();
            hits.append(inter);

            point pt = inter.point() + smallVec[pointI];

            while (((pt-start[pointI])&dirVec[pointI]) <= magSqrDirVec[pointI])
            {
                // See if any intersection between pt and end
                pointIndexHit inter = findLine(pt, end[pointI]);

                // Check for not hit or hit same face as before (can happen
                // if vector along surface of face)
                if
                (
                    !inter.hit()
                 || (inter.index() == hits.last().index())
                )
                {
                    break;
                }
                hits.append(inter);

                pt = inter.point() + smallVec[pointI];
            }

            info[pointI].transfer(hits);
        }
        else
        {
            info[pointI].clear();
        }
    }
}


void Foam::searchableRotatedBox::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableRotatedBox::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    // searchableBox does not use hitPoints so no need to transform
    box_.getNormal(info, normal);

    normal = transform_.globalVector(normal);
}


void Foam::searchableRotatedBox::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    box_.getVolumeType(transform_.localPosition(points), volType);
}


// ************************************************************************* //
