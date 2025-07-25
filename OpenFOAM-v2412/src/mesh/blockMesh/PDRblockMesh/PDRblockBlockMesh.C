/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "PDRblock.H"
#include "ListOps.H"
#include "gradingDescriptors.H"
#include "objectRegistry.H"
#include "Time.H"
#include "IOdictionary.H"
#include "Fstream.H"
#include "OTstream.H"
#include "edgeHashes.H"

#include "cellModel.H"
#include "blockMesh.H"
#include "polyMesh.H"
#include "searchableSphere.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Name for the projection geometry
static const word projKeyword("project");

// Name for the projection geometry
static const word projGeomName("sphere");


//- Calculate geometric ratio from relative ratio
inline scalar relativeToGeometricRatio
(
    const scalar expRatio,
    const label nDiv
)
{
    return nDiv > 1 ? pow(expRatio, (nDiv - 1)) : 1.0;
}


// Output space-separated flat list. No size prefix.
template<class T>
static Ostream& outputFlatList(Ostream& os, const UList<T>& list)
{
    os  << token::BEGIN_LIST;
    label i = 0;
    for (const label val : list)
    {
        if (i++) os << token::SPACE;
        os  << val;
    }
    os  << token::END_LIST;

    return os;
}


// Begin indent list
static inline Ostream& begIndentList(Ostream& os)
{
    os  << indent << incrIndent << token::BEGIN_LIST << nl;
    return os;
}

// End indent list
static inline Ostream& endIndentList(Ostream& os)
{
    os  << decrIndent << indent << token::END_LIST;
    return os;
}


// Output list contents (newline separated) indented.
template<class T>
static Ostream& outputIndent(Ostream& os, const UList<T>& list)
{
    for (const T& val : list)
    {
        os  << indent << val << nl;
    }
    return os;
}


//
static Ostream& serializeHex
(
    Ostream& os,
    const labelUList& hexVerts,
    const labelVector& hexCount,
    const Vector<gradingDescriptors> hexGrade,
    const word& zoneName = word::null
)
{
    os  << indent << cellModel::modelNames[cellModel::HEX] << token::SPACE;
    outputFlatList(os, hexVerts);

    if (!zoneName.empty())
    {
        os  << token::SPACE << zoneName;
    }

    os  << token::SPACE << hexCount << nl
        << indent << word("edgeGrading") << nl;

    begIndentList(os);

    // Grading (x/y/z)
    for (const gradingDescriptors& gds : hexGrade)
    {
        begIndentList(os);
        outputIndent(os, gds);
        endIndentList(os) << nl;
    }

    endIndentList(os) << nl;
    return os;
}


// Generate list with entries:
//
//    project (x y z) (geometry)
//
static Ostream& serializeProjectPoints
(
    Ostream& os,
    const UList<point>& list
)
{
    for (const point& p : list)
    {
        os  << indent << projKeyword << token::SPACE
            << p
            << token::SPACE
            << token::BEGIN_LIST << projGeomName << token::END_LIST << nl;
    }

    return os;
}


// Generate entry:
//
//    project (beg end) (geometry)
//
static Ostream& serializeProjectEdge
(
    Ostream& os,
    const edge& e
)
{
    os  << indent << projKeyword << token::SPACE;

    if (e.is_sorted())
    {
        os  << e.first() << token::SPACE << e.second();
    }
    else
    {
        os  << e.second() << token::SPACE << e.first();
    }

    os  << token::SPACE
        << token::BEGIN_LIST << projGeomName << token::END_LIST << nl;

    return os;
}


// Generate entry:
//
//    (0 1 2 ..)
//
static Ostream& serializeFace
(
    Ostream& os,
    const face& list
)
{
    os  << indent;
    outputFlatList(os, list);
    os  << nl;
    return os;
}


// Generate entry:
//
//    project (0 1 2 ..) geometry
//
static Ostream& serializeProjectFace
(
    Ostream& os,
    const face& list
)
{
    os  << indent << projKeyword << token::SPACE;
    outputFlatList(os, list);
    os  << token::SPACE << projGeomName << nl;

    return os;
}

} // namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::PDRblock::blockMeshDict
(
    Ostream& os,
    const bool withHeader
) const
{
    if (withHeader)
    {
        // Use dummy time for fake objectRegistry
        autoPtr<Time> dummyTimePtr(Time::New());

        IOdictionary iodict
        (
            IOobject
            (
                "blockMeshDict",
                dummyTimePtr->system(),
                *dummyTimePtr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        iodict.writeHeader(os);
    }

    const cellModel& hex = cellModel::ref(cellModel::HEX);

    // The mesh topology will normally be an O-grid with a central (inner)
    // block and 6 outer blocks.
    // The inner block is described by (0 1 2 3 4 5 6 7).
    // The additional points for the outer region: (8 9 19 11 12 13 14 15)

    // The outer blocks will be addressed according to their
    // placement w.r.t. the inner block, and defined such that they retain
    // the same global orientation as the inner block.
    // For example, face 0 of all blocks will be on the logical x-min side
    // of the mesh.

    // List of hex vertices, ordered with outer blocks first to allow
    // direct addressing by their logical position.
    const FixedList<labelList, 8>
        hexVerts
        ({
            {8, 0,  3, 11, 12,  4,  7, 15},     // x-min block (face 0)
            {1, 9, 10,  2,  5, 13, 14,  6},     // x-max block (face 1)
            {8, 9,  1,  0, 12, 13,  5,  4},     // y-min block (face 2)
            {3, 2, 10, 11,  7,  6, 14, 15},     // y-max block (face 3)
            {8, 9, 10, 11,  0,  1,  2,  3},     // z-min block (face 4)
            {4, 5,  6,  7, 12, 13, 14, 15},     // z-max block (face 5)
            {0, 1,  2,  3,  4,  5,  6,  7},     // Inner box description
            {8, 9, 10, 11, 12, 13, 14, 15},     // Outer box description
        });

    // The face or the logical block index (for the O-grid)
    enum faceIndex
    {
        X_Min = 0,
        X_Max = 1,
        Y_Min = 2,
        Y_Max = 3,
        Z_Min = 4,
        Z_Max = 5,
        Inner_Block = 6,    // Inner block description
        Outer_Block = 7,    // Outer bounding box description
    };

    // Lists of block/face for outside and ground faces
    DynamicList<labelPair> outerFaces(8);
    DynamicList<labelPair> groundFaces(8);

    // We handle a few fixed topology configurations

    enum blockTopologyType
    {
        INNER_ONLY = 0,   // No outer region
        EXTENDED   = 1,   // Outer created by extending inner region
        CLIP_BOTTOM = 5,  // Outer O-grid on 5 sides
        FULL_OUTER = 6    // Outer O-grid on 6 sides
    };



    // Expansion ratios need conversion from relative to geometric
    const bool useRelToGeom =
        (expansionType::EXPAND_RATIO == outer_.expandType_);


    // Physical dimensions

    Vector<gridControl> ctrl(control_);
    boundBox innerCorners(bounds(ctrl.x(), ctrl.y(), ctrl.z()));

    boundBox outerCorners;


    point radialCentre(innerCorners.centre());
    vector radialSizes(0.5*innerCorners.span());

    blockTopologyType outerTopology = INNER_ONLY;


    if (outer_.active())
    {
        outerTopology = FULL_OUTER;

        // Convert from relative size
        radialSizes.x() *= outer_.relSize_.x();
        radialSizes.y() *= outer_.relSize_.y();
        radialSizes.z() *= min(outer_.relSize_.x(), outer_.relSize_.y());

        if (outer_.onGround())
        {
            outerTopology = CLIP_BOTTOM;
            radialCentre.z() = innerCorners.min().z();
            radialSizes.z() *= 2;
        }

        // Corners for a box
        outerCorners.min() = radialCentre - radialSizes;
        outerCorners.max() = radialCentre + radialSizes;

        if (outer_.onGround())
        {
            outerCorners.min().z() = innerCorners.min().z();
        }


        if (outer_.isSphere())
        {
            // For spheroid projection, don't trust that blockMesh does it
            // properly. Give some reasonable estimates of the corners.

            // Use dummy Time for objectRegistry
            autoPtr<Time> dummyTimePtr(Time::New());

            const IOobject io
            (
                "sphere",
                *dummyTimePtr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            searchableSphere sphere(io, radialCentre, radialSizes);

            pointField queries(2);
            queries[0] = outerCorners.min();
            queries[1] = outerCorners.max();

            List<pointIndexHit> hits;
            sphere.findNearest
            (
                queries,
                scalarField(2, GREAT),
                hits
            );

            outerCorners.min() = hits[0].hitPoint();
            outerCorners.max() = hits[1].hitPoint();
        }
        else if (outerControl::OUTER_EXTEND == outer_.type_)
        {
            outerTopology = EXTENDED;

            // Extend the inner block
            label outerCount;
            scalar expRatio;

            outerCount = outer_.nCells_.x();
            expRatio = outer_.expansion_.x();
            if (useRelToGeom)
            {
                expRatio = relativeToGeometricRatio(expRatio, outerCount);
            }

            ctrl.x().prepend(outerCorners.min().x(), outerCount, -expRatio);
            ctrl.x().append(outerCorners.max().x(), outerCount, expRatio);


            outerCount = outer_.nCells_.y();
            expRatio = outer_.expansion_.y();
            if (useRelToGeom)
            {
                expRatio = relativeToGeometricRatio(expRatio, outerCount);
            }

            ctrl.y().prepend(outerCorners.min().y(), outerCount, -expRatio);
            ctrl.y().append(outerCorners.max().y(), outerCount, expRatio);

            outerCount = max(outer_.nCells_.x(), outer_.nCells_.y());
            expRatio = min(outer_.expansion_.x(), outer_.expansion_.y());
            if (useRelToGeom)
            {
                expRatio = relativeToGeometricRatio(expRatio, outerCount);
            }

            if (!outer_.onGround())
            {
                ctrl.z().prepend(outerCorners.min().z(), outerCount, -expRatio);
            }
            ctrl.z().append(outerCorners.max().z(), outerCount, expRatio);

            // Update corners
            innerCorners = bounds(ctrl.x(), ctrl.y(), ctrl.z());
            outerCorners = innerCorners;
        }
    }


    const Vector<gradingDescriptors> innerGrading(grading(ctrl));
    const labelVector innerCount(sizes(ctrl));

    labelVector hexCount;
    Vector<gradingDescriptors> hexGrade;


    const label radialCount = outer_.nCells_.x();
    scalar expRatio = outer_.expansion_.x();

    if (useRelToGeom)
    {
        expRatio = relativeToGeometricRatio(expRatio, radialCount);
    }

    const gradingDescriptors radialInward
    (
        gradingDescriptor{-expRatio}
    );

    const gradingDescriptors radialOutward
    (
        gradingDescriptor{expRatio}
    );


    if (EXTENDED == outerTopology)
    {
        // The inner block is extended to become the outer faces
        outerFaces.append
        ({
            labelPair(Inner_Block, X_Min),
            labelPair(Inner_Block, X_Max),
            labelPair(Inner_Block, Y_Min),
            labelPair(Inner_Block, Y_Max),
            labelPair(Inner_Block, Z_Max)
        });

        // The ground faces vs outside faces
        if (outer_.onGround())
        {
            groundFaces.append
            (
                labelPair(Inner_Block, Z_Min)
            );
        }
        else
        {
            outerFaces.append
            (
                labelPair(Inner_Block, Z_Min)
            );
        }
    }
    else if (CLIP_BOTTOM == outerTopology || FULL_OUTER == outerTopology)
    {
        // The outside faces
        outerFaces.append
        ({
            labelPair(X_Min, X_Min),
            labelPair(X_Max, X_Max),
            labelPair(Y_Min, Y_Min),
            labelPair(Y_Max, Y_Max),
            labelPair(Z_Max, Z_Max)
        });

        // The ground faces
        if (CLIP_BOTTOM == outerTopology)
        {
            groundFaces.append
            ({
                labelPair(X_Min, Z_Min),
                labelPair(X_Max, Z_Min),
                labelPair(Y_Min, Z_Min),
                labelPair(Y_Max, Z_Min),
                // Note: {Z_Min, Z_Min} will not exist
                labelPair(Inner_Block, Z_Min)
            });
        }
        else
        {
            outerFaces.append
            (
                labelPair(Z_Min, Z_Min)
            );
        }
    }


    if (outer_.isSphere())
    {
        os.beginBlock("geometry");
        {
            os.beginBlock(projGeomName);
            {
                os.writeEntry("type", "sphere");
                os.writeEntry("origin", radialCentre);
                os.writeEntry("radius", radialSizes);
            }
            os.endBlock();
        }
        os.endBlock();
    }

    // vertices
    {
        os  << nl << word("vertices") << nl;
        begIndentList(os);

        pointField corners(innerCorners.points());

        // inner
        outputIndent(os, corners);

        // outer
        if (CLIP_BOTTOM == outerTopology || FULL_OUTER == outerTopology)
        {
            corners = outerCorners.points();

            if (outer_.isSphere())
            {
                serializeProjectPoints(os, corners);
            }
            else
            {
                outputIndent(os, corners);
            }
        }

        endIndentList(os);
        os.endEntry();
    }


    // blocks
    {
        word innerZoneName = "inner";
        if (INNER_ONLY == outerTopology || EXTENDED == outerTopology)
        {
            innerZoneName.clear();
        }

        os  << nl << word("blocks") << nl;
        begIndentList(os);

        // Inner block
        hexCount = innerCount;

        serializeHex
        (
            os,
            hexVerts[Inner_Block],
            hexCount,
            innerGrading,
            innerZoneName
        );

        // outer
        if (CLIP_BOTTOM == outerTopology || FULL_OUTER == outerTopology)
        {
            // Radial direction = X
            hexCount = innerCount;
            hexGrade = innerGrading;

            hexCount.x() = radialCount;

            // Face 0: x-min
            {
                hexGrade.x() = radialInward;
                serializeHex(os, hexVerts[X_Min], hexCount, hexGrade);
            }
            // Face 1: x-max
            {
                hexGrade.x() = radialOutward;
                serializeHex(os, hexVerts[X_Max], hexCount, hexGrade);
            }


            // Radial direction = Y
            hexCount = innerCount;
            hexGrade = innerGrading;

            hexCount.y() = radialCount;

            // Face 2: y-min
            {
                hexGrade.y() = radialInward;
                serializeHex(os, hexVerts[Y_Min], hexCount, hexGrade);
            }
            // Face 3: y-max
            {
                hexGrade.y() = radialOutward;
                serializeHex(os, hexVerts[Y_Max], hexCount, hexGrade);
            }


            // Radial direction = Z
            hexCount = innerCount;
            hexGrade = innerGrading;

            hexCount.z() = radialCount;

            // Face 4: z-min
            if (!outer_.onGround())
            {
                hexGrade.z() = radialInward;
                serializeHex(os, hexVerts[Z_Min], hexCount, hexGrade);
            }
            // Face 5: z-max
            {
                hexGrade.z() = radialOutward;
                serializeHex(os, hexVerts[Z_Max], hexCount, hexGrade);
            }
        }

        endIndentList(os);
        os.endEntry();
    }


    // edges
    {
        os  << nl << word("edges") << nl;
        begIndentList(os);

        if (outer_.isSphere() && outerFaces.size())
        {
            // Edges for the outer face of the block
            edgeHashSet projEdges(32);

            for (const labelPair& pr : outerFaces)
            {
                projEdges.insert
                (
                    hex.face(pr.second(), hexVerts[pr.first()]).edges()
                );
            }

            for (const edge& e : projEdges.sortedToc())
            {
                serializeProjectEdge(os, e);
            }
        }

        endIndentList(os);
        os.endEntry();
    }


    // faces
    {
        os  << nl << word("faces") << nl;
        begIndentList(os);

        if (outer_.isSphere() && outerFaces.size())
        {
            for (const labelPair& pr : outerFaces)
            {
                serializeProjectFace
                (
                    os,
                    hex.face(pr.second(), hexVerts[pr.first()])
                );
            }
        }

        endIndentList(os);
        os.endEntry();
    }


    // boundary
    {
        os  << nl << word("boundary") << nl;
        begIndentList(os);

        // outer
        {
            os.beginBlock("outer");
            os.writeEntry("type", word("patch"));

            os  << indent << word("faces") << nl;
            begIndentList(os);

            for (const labelPair& pr : outerFaces)
            {
                serializeFace
                (
                    os,
                    hex.face(pr.second(), hexVerts[pr.first()])
                );
            }

            endIndentList(os);
            os.endEntry();

            os.endBlock();
        }

        if (outer_.onGround())
        {
            os.beginBlock("ground");
            os.writeEntry("type", word("wall"));

            os  << indent << word("faces") << nl;
            begIndentList(os);

            for (const labelPair& pr : groundFaces)
            {
                serializeFace
                (
                    os,
                    hex.face(pr.second(), hexVerts[pr.first()])
                );
            }

            endIndentList(os);
            os.endEntry();

            os.endBlock();
        }

        endIndentList(os);
        os.endEntry();
    }


    if (withHeader)
    {
        IOobject::writeEndDivider(os);
    }

    return os;
}


Foam::dictionary Foam::PDRblock::blockMeshDict() const
{
    OTstream os;
    blockMeshDict(os);

    ITstream is;
    is.transfer(os.tokens());

    return dictionary(is);
}


void Foam::PDRblock::writeBlockMeshDict(const IOobject& io) const
{
    IOdictionary iodict
    (
        IOobject
        (
            io.name(),
            io.db().time().system(),
            io.local(),
            io.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    OFstream os(iodict.objectPath());

    Info<< nl
        << "Generate blockMeshDict: "
        << iodict.db().time().relativePath(os.name()) << endl;

    // More precision (for points data)
    os.precision(Foam::max(10u, IOstream::defaultPrecision()));

    iodict.writeHeader(os);

    // Just like writeData, but without copying beforehand
    this->blockMeshDict(os);

    IOobject::writeEndDivider(os);
}


Foam::autoPtr<Foam::blockMesh>
Foam::PDRblock::createBlockMesh(const IOobject& io) const
{
    IOdictionary iodict
    (
        IOobject
        (
            "blockMeshDict.PDRblockMesh",
            io.db().time().system(),
            io.local(),
            io.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        blockMeshDict()
    );

    return autoPtr<blockMesh>::New(iodict);
}


Foam::autoPtr<Foam::polyMesh>
Foam::PDRblock::meshBlockMesh(const IOobject& io) const
{
    const bool oldVerbose = blockMesh::verboseOutput;
    blockMesh::verboseOutput = false;

    autoPtr<polyMesh> meshPtr(createBlockMesh(io)->mesh(io));

    blockMesh::verboseOutput = oldVerbose;

    // This is a bit ugly.
    // For extend, we still wish to have an 'inner' cellZone,
    // but we meshed the entirety.

    if
    (
        outerControl::OUTER_EXTEND == outer_.type_
     && meshPtr->cellZones().empty()
    )
    {
        const boundBox innerBox
        (
            bounds(control_.x(), control_.y(), control_.z())
        );

        const label nZoneCellsMax =
        (
            control_.x().nCells()
          * control_.y().nCells()
          * control_.z().nCells()
        );


        polyMesh& pmesh = *meshPtr;

        List<cellZone*> cz(1);
        cz[0] = new cellZone
        (
            "inner",
            labelList(nZoneCellsMax),
            0, // zonei
            pmesh.cellZones()
        );

        cellZone& innerZone = *(cz[0]);

        const vectorField& cc = pmesh.cellCentres();

        label nZoneCells = 0;

        for
        (
            label celli = 0;
            celli < cc.size() && nZoneCells < nZoneCellsMax;
            ++celli
        )
        {
            if (innerBox.contains(cc[celli]))
            {
                innerZone[nZoneCells] = celli;
                ++nZoneCells;
            }
        }

        innerZone.resize(nZoneCells);

        pmesh.pointZones().clear();
        pmesh.faceZones().clear();
        pmesh.cellZones().clear();
        pmesh.addZones(List<pointZone*>(), List<faceZone*>(), cz);
    }

    return meshPtr;
}


// ************************************************************************* //
