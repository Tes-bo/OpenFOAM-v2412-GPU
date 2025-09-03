/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "MappedFileFilterField.H"
#include "oneField.H"
#include "MeshedSurfaces.H"
#include "MinMax.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{
    defineTypeNameAndDebug(FilterField, 0);
}
}


const Foam::Enum
<
    Foam::PatchFunction1Types::FilterField::RBF_type
>
Foam::PatchFunction1Types::FilterField::RBF_typeNames_
{
    { RBF_type::RBF_linear, "linear" },
    { RBF_type::RBF_quadratic, "quadratic" },
    { RBF_type::RBF_linear, "default" },
};


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Linear basis function: 1 - (r/r0)
static inline scalar linearWeight
(
    const point& p,
    const point& p0,
    const scalar radiusSqr
)
{
    return (1 - ::sqrt(p.distSqr(p0) / radiusSqr));
}


//- Quadratic basis function: 1 - (r/r0)^2
static inline scalar quadraticWeight
(
    const point& p,
    const point& p0,
    const scalar radiusSqr
)
{
    return (1 - p.distSqr(p0) / radiusSqr);
}


//- Construct search tree for points
static autoPtr<indexedOctree<treeDataPoint>>
createTree(const pointField& points)
{
    // Avoid bounding box error in case of 2D cases
    treeBoundBox bb(points);
    bb.grow(1e-4);

    const int debug = PatchFunction1Types::FilterField::debug;

    const int oldDebug = indexedOctreeBase::debug;

    if (debug & 2)
    {
        indexedOctreeBase::debug = 1;
    }

    auto treePtr = autoPtr<indexedOctree<treeDataPoint>>::New
    (
        treeDataPoint(points),
        bb,             // overall search domain
        points.size(),  // maxLevel
        16,             // leafsize
        1.0             // duplicity
    );

    indexedOctreeBase::debug = oldDebug;

    if (debug & 2)
    {
        const auto& tree = treePtr();

        OFstream os("indexedOctree.obj");
        tree.writeOBJ(os);
        Info<< "=================" << endl;
        Info<< "have " << tree.nodes().size() << " nodes" << nl;
        Info<< "=================" << endl;
    }

    return treePtr;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class TreeType,
    class RadiusSqrOp,
    class BasisFunctionOp,
    class PointWeightFieldType
>
void Foam::PatchFunction1Types::FilterField::buildWeightsImpl
(
    const TreeType& tree,
    const RadiusSqrOp& radiusSqrOp,
    const BasisFunctionOp& basisFuncOp,
    const PointWeightFieldType& pointWeightFld
)
{
    tmp<pointField> tpoints = tree.shapes().centres();

    const auto& points = tpoints();

    const label nPoints = points.size();

    weights_.clear();
    addressing_.clear();

    weights_.resize(nPoints);
    addressing_.resize(nPoints);

    labelHashSet neighbours(2*128);

    // Catch degenerate case where weight/addressing not actually needed
    bool usesNeighbours = false;

    for (label pointi = 0; pointi < nPoints; ++pointi)
    {
        const point& p0 = points[pointi];
        auto& currAddr = addressing_[pointi];
        auto& currWeight = weights_[pointi];

        // Search with local sphere
        const scalar radiusSqr = radiusSqrOp(pointi);

        tree.findSphere(p0, radiusSqr, neighbours);

        // Paranoid, enforce identity weighting as a minimum
        if (neighbours.size() < 2)
        {
            currAddr.resize(1, pointi);
            currWeight.resize(1, scalar(1));
            continue;
        }

        usesNeighbours = true;

        currAddr = neighbours.sortedToc();
        currWeight.resize_nocopy(currAddr.size());

        scalar sumWeight = 0;

        forAll(currAddr, i)
        {
            const point& p = points[currAddr[i]];

            currWeight[i] = basisFuncOp(p, p0, radiusSqr);

            // Eg, face area weighting.
            // - cast required for oneField()
            currWeight[i] *= static_cast<scalar>(pointWeightFld[i]);
            sumWeight += currWeight[i];
        }

        for (auto& w : currWeight)
        {
            w /= sumWeight;
        }
    }

    if (!usesNeighbours)
    {
        // No addressing/weighting required
        reset();
    }

    if (debug && !addressing_.empty())
    {
        labelMinMax limits(addressing_[0].size());

        label total = 0;

        for (const labelList& addr : addressing_)
        {
            const label n = addr.size();

            limits.add(n);
            total += n;
        }

        Pout<< "Weight neighbours: min=" << limits.min()
            << " avg=" << (total / addressing_.size())
            << " max=" << limits.max() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PatchFunction1Types::FilterField::FilterField
(
    const pointField& points,
    const scalar radius,
    const RBF_type interp
)
{
    reset(points, radius, interp);
}


Foam::PatchFunction1Types::FilterField::FilterField
(
    const meshedSurface& geom,
    const scalar radius,
    const bool relative,
    const RBF_type interp
)
{
    reset(geom, radius, relative, interp);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PatchFunction1Types::FilterField::reset()
{
    addressing_.clear();
    weights_.clear();
}


void Foam::PatchFunction1Types::FilterField::reset
(
    const pointField& points,
    const scalar radius,
    const RBF_type interp
)
{
    scalarField dummyWeights;

    if (debug)
    {
        Pout<< "Apply " << RBF_typeNames_[interp] << " filter,"
            << " radius=" << radius << nl
            << "Create tree..." << endl;
    }

    autoPtr<indexedOctree<treeDataPoint>> treePtr
        = createTree(points);

    const scalar radiusSqr = sqr(radius);

    auto weightFunc = linearWeight;

    if (interp == RBF_type::RBF_quadratic)
    {
        weightFunc = quadraticWeight;
    }


    // Use (constant) absolute size
    buildWeightsImpl
    (
        treePtr(),
        [=](label index) -> scalar { return radiusSqr; },
        weightFunc,
        oneField()
    );
}


void Foam::PatchFunction1Types::FilterField::reset
(
    const meshedSurface& geom,
    const scalar radius,
    const bool relative,
    const RBF_type interp
)
{
    if (!geom.nFaces())
    {
        // Received geometry without any faces eg, boundaryData
        reset(geom.points(), radius, interp);
        return;
    }

    if (debug)
    {
        Pout<< "Apply " << RBF_typeNames_[interp] << " filter,"
            << (relative ? " relative" : "") << " radius=" << radius << nl
            << "Create tree..." << endl;
    }

    autoPtr<indexedOctree<treeDataPoint>> treePtr
        = createTree(geom.faceCentres());

    const scalar radiusSqr = sqr(radius);

    auto weightFunc = linearWeight;

    if (interp == RBF_type::RBF_quadratic)
    {
        weightFunc = quadraticWeight;
    }


    if (relative)
    {
        // Use relative size
        buildWeightsImpl
        (
            treePtr(),
            [&](label index) -> scalar
            {
                return (radiusSqr * geom.sphere(index));
            },
            weightFunc,
            geom.magSf()
        );
    }
    else
    {
        // Use (constant) absolute size
        buildWeightsImpl
        (
            treePtr(),
            [=](label index) -> scalar { return radiusSqr; },
            weightFunc,
            geom.magSf()
        );
    }
}


// ************************************************************************* //
