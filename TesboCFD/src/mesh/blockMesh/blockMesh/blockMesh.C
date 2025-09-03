/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "blockMesh.H"
#include "transform.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitch(blockMesh, 0);
}

bool Foam::blockMesh::verboseOutput = true;


const Foam::Enum<Foam::blockMesh::mergeStrategy>
Foam::blockMesh::strategyNames_
({
    { mergeStrategy::MERGE_TOPOLOGY, "topology" },
    { mergeStrategy::MERGE_POINTS, "points" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Define/refine the verbosity level
// Command-line options have precedence over dictionary setting
static int getVerbosity(const dictionary& dict, int verbosity)
{
    if (verbosity < 0)
    {
        // Forced as 'off'
        verbosity = 0;
    }
    else if (!verbosity)
    {
        // Not specified: use dictionary value or static default
        verbosity = dict.getOrDefault("verbose", blockMesh::verboseOutput);
    }

    return verbosity;
}


// Process dictionary entry with single scalar or vector quantity
// - return 0 if scaling is not needed. Eg, Not found or is unity
// - return 1 for uniform scaling
// - return 3 for non-uniform scaling

int readScaling(const entry *eptr, vector& scale)
{
    int nCmpt = 0;

    if (!eptr)
    {
        return nCmpt;
    }

    ITstream& is = eptr->stream();

    if (is.peek().isNumber())
    {
        // scalar value
        scalar val;
        is >> val;

        if ((val > 0) && !equal(val, 1))
        {
            // Uniform scaling
            nCmpt = 1;
            scale = vector::uniform(val);
        }
    }
    else
    {
        // vector value
        is >> scale;

        bool nonUnity = false;
        for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
        {
            if (scale[cmpt] <= 0)
            {
                scale[cmpt] = 1;
            }
            else if (!equal(scale[cmpt], 1))
            {
                nonUnity = true;
            }
        }

        if (nonUnity)
        {
            if (equal(scale.x(), scale.y()) && equal(scale.x(), scale.z()))
            {
                // Uniform scaling
                nCmpt = 1;
            }
            else
            {
                nCmpt = 3;
            }
        }
    }

    eptr->checkITstream(is);

    return nCmpt;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::blockMesh::readPointTransforms(const dictionary& dict)
{
    transformType_ = transformTypes::NO_TRANSFORM;

    // Optional (cartesian) coordinate system transform
    auto csysPtr = coordinateSystem::NewIfPresent(dict, "transform");

    if (csysPtr)
    {
        transform_ = csysPtr();

        // Non-zero origin?
        if (magSqr(transform_.origin()) > ROOTVSMALL)
        {
            transformType_ |= transformTypes::TRANSLATION;
        }

        // Non-identity rotation?
        if (!transform_.R().is_identity(ROOTVSMALL))
        {
            transformType_ |= transformTypes::ROTATION;
        }
    }
    else
    {
        transform_.clear();
    }


    // Optional 'prescale' factor.
    {
        prescaling_ = vector::uniform(1);

        const int scaleType = readScaling
        (
            dict.findEntry("prescale", keyType::LITERAL),
            prescaling_
        );

        if (scaleType == 1)
        {
            transformType_ |= transformTypes::PRESCALING;
        }
        else if (scaleType == 3)
        {
            transformType_ |= transformTypes::PRESCALING3;
        }
    }


    // Optional 'scale' factor. Was 'convertToMeters' until OCT-2008
    {
        scaling_ = vector::uniform(1);

        const int scaleType = readScaling
        (
            // Mark as changed from 2010 onwards
            dict.findCompat
            (
                "scale",
                {{"convertToMeters", 1012}},
                keyType::LITERAL
            ),
            scaling_
        );

        if (scaleType == 1)
        {
            transformType_ |= transformTypes::SCALING;
        }
        else if (scaleType == 3)
        {
            transformType_ |= transformTypes::SCALING3;
        }
    }

    return bool(transformType_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMesh::blockMesh
(
    const IOdictionary& dict,
    const word& regionName,
    mergeStrategy strategy,
    int verbosity
)
:
    meshDict_(dict),
    verbose_(getVerbosity(dict, verbosity)),
    checkFaceCorrespondence_
    (
        meshDict_.getOrDefault("checkFaceCorrespondence", true)
    ),
    mergeStrategy_(strategy),
    transformType_(transformTypes::NO_TRANSFORM),
    geometry_
    (
        IOobject
        (
            "geometry",                 // dummy name
            meshDict_.time().constant(), // instance
            "geometry",                 // local
            meshDict_.time(),           // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshDict_.found("geometry")
      ? meshDict_.subDict("geometry")
      : dictionary(),
        true
    ),
    blockVertices_
    (
        meshDict_.lookup("vertices"),
        blockVertex::iNew(meshDict_, geometry_)
    ),
    vertices_(Foam::vertices(blockVertices_)),
    prescaling_(vector::uniform(1)),
    scaling_(vector::uniform(1)),
    transform_(),
    topologyPtr_(createTopology(meshDict_, regionName))
{
    if (mergeStrategy_ == mergeStrategy::DEFAULT_MERGE)
    {
        strategyNames_.readIfPresent("mergeType", meshDict_, mergeStrategy_);

        // Warn about fairly obscure old "fastMerge" option?

        if
        (
            mergeStrategy_ == mergeStrategy::DEFAULT_MERGE
         && checkDegenerate()
        )
        {
            Info<< nl
                << "Detected collapsed blocks "
                << "- using merge points instead of merge topology" << nl
                << endl;

            mergeStrategy_ = mergeStrategy::MERGE_POINTS;
        }
    }

    if (mergeStrategy_ == mergeStrategy::MERGE_POINTS)
    {
        // MERGE_POINTS
        calcGeometricalMerge();
    }
    else
    {
        // MERGE_TOPOLOGY
        calcTopologicalMerge();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::blockMesh::good() const noexcept
{
    return bool(topologyPtr_);
}


int Foam::blockMesh::verbose() const noexcept
{
    return verbose_;
}


int Foam::blockMesh::verbose(const int level) noexcept
{
    int old(verbose_);
    verbose_ = level;
    return old;
}


const Foam::pointField& Foam::blockMesh::vertices() const noexcept
{
    return vertices_;
}


Foam::tmp<Foam::pointField>
Foam::blockMesh::vertices(bool applyTransform) const
{
    if (applyTransform && hasPointTransforms())
    {
        auto tpts = tmp<pointField>::New(vertices_);

        inplacePointTransforms(tpts.ref());

        return tpts;
    }

    return vertices_;
}


Foam::PtrList<Foam::dictionary> Foam::blockMesh::patchDicts() const
{
    const polyPatchList& topoPatches = topology().boundaryMesh();

    PtrList<dictionary> patchDicts(topoPatches.size());

    OCharStream os;
    ISpanStream is;

    forAll(topoPatches, patchi)
    {
        os.rewind();
        topoPatches[patchi].write(os);

        is.reset(os.view());
        patchDicts.emplace_set(patchi, is);
    }
    return patchDicts;
}


const Foam::pointField& Foam::blockMesh::points() const
{
    if (points_.empty())
    {
        createPoints();
    }

    return points_;
}


const Foam::cellShapeList& Foam::blockMesh::cells() const
{
    if (cells_.empty())
    {
        createCells();
    }

    return cells_;
}


const Foam::faceListList& Foam::blockMesh::patches() const
{
    if (patches_.empty())
    {
        createPatches();
    }

    return patches_;
}


Foam::wordList Foam::blockMesh::patchNames() const
{
    return topology().boundaryMesh().names();
}


//Foam::wordList Foam::blockMesh::patchTypes() const
//{
//    return topology().boundaryMesh().types();
//}
//
//
//Foam::wordList Foam::blockMesh::patchPhysicalTypes() const
//{
//    return topology().boundaryMesh().physicalTypes();
//}


Foam::label Foam::blockMesh::numZonedBlocks() const
{
    const blockList& blocks = *this;

    label count = 0;

    for (const block& blk : blocks)
    {
        if (blk.zoneName().size())
        {
            ++count;
        }
    }

    return count;
}


bool Foam::blockMesh::hasPointTransforms() const noexcept
{
    return bool(transformType_);
}


bool Foam::blockMesh::inplacePointTransforms(pointField& pts) const
{
    if (!transformType_)
    {
        return false;
    }

    if (transformType_ & transformTypes::PRESCALING)
    {
        for (point& p : pts)
        {
            p *= prescaling_.x();
        }
    }
    else if (transformType_ & transformTypes::PRESCALING3)
    {
        for (point& p : pts)
        {
            p = cmptMultiply(p, prescaling_);
        }
    }

    if (transformType_ & transformTypes::ROTATION)
    {
        const tensor rot(transform_.R());

        if (transformType_ & transformTypes::TRANSLATION)
        {
            const point origin(transform_.origin());

            for (point& p : pts)
            {
                p = Foam::transform(rot, p) + origin;
            }
        }
        else
        {
            for (point& p : pts)
            {
                p = Foam::transform(rot, p);
            }
        }
    }
    else if (transformType_ & transformTypes::TRANSLATION)
    {
        const point origin(transform_.origin());

        for (point& p : pts)
        {
            p += origin;
        }
    }

    if (transformType_ & transformTypes::SCALING)
    {
        for (point& p : pts)
        {
            p *= scaling_.x();
        }
    }
    else if (transformType_ & transformTypes::SCALING3)
    {
        for (point& p : pts)
        {
            p = cmptMultiply(p, scaling_);
        }
    }

    return true;
}


Foam::tmp<Foam::pointField>
Foam::blockMesh::globalPosition(const pointField& localPoints) const
{
    if (hasPointTransforms())
    {
        auto tpts = tmp<pointField>::New(localPoints);

        inplacePointTransforms(tpts.ref());

        return tpts;
    }
    else
    {
        return localPoints;
    }
}


// ************************************************************************* //
