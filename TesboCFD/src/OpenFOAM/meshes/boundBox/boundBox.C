/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "boundBox.H"
#include "PstreamReduceOps.H"
#include "plane.H"
#include "hexCell.H"
#include "triangle.H"
#include "MinMax.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::boundBox Foam::boundBox::greatBox
(
    point::uniform(-ROOTVGREAT),
    point::uniform(ROOTVGREAT)
);

const Foam::boundBox Foam::boundBox::invertedBox
(
    point::uniform(ROOTVGREAT),
    point::uniform(-ROOTVGREAT)
);

const Foam::FixedList<Foam::vector, 6> Foam::boundBox::faceNormals
({
    vector(-1,  0,  0), // 0: x-min, left
    vector( 1,  0,  0), // 1: x-max, right
    vector( 0, -1,  0), // 2: y-min, bottom
    vector( 0,  1,  0), // 3: y-max, top
    vector( 0,  0, -1), // 4: z-min, back
    vector( 0,  0,  1)  // 5: z-max, front
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::faceList& Foam::boundBox::hexFaces()
{
    return hexCell::modelFaces();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const boundBox& bb, const bool doReduce)
:
    boundBox(bb)
{
    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox(const UList<point>& points, bool doReduce)
:
    boundBox()
{
    add(points);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox(const tmp<pointField>& tpoints, bool doReduce)
:
    boundBox()
{
    add(tpoints);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox
(
    const UList<point>& points,
    const labelUList& indices,
    bool doReduce
)
:
    boundBox()
{
    add(points, indices);

    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundBox::hexCorners() const
{
    auto tpts = tmp<pointField>::New(boundBox::nPoints());
    auto& pts = tpts.ref();

    pts[0] = hexCorner<0>();
    pts[1] = hexCorner<1>();
    pts[2] = hexCorner<2>();
    pts[3] = hexCorner<3>();
    pts[4] = hexCorner<4>();
    pts[5] = hexCorner<5>();
    pts[6] = hexCorner<6>();
    pts[7] = hexCorner<7>();

    return tpts;
}


Foam::tmp<Foam::pointField> Foam::boundBox::faceCentres() const
{
    auto tpts = tmp<pointField>::New(boundBox::nFaces());
    auto& pts = tpts.ref();

    for (direction facei = 0; facei < boundBox::nFaces(); ++facei)
    {
        pts[facei] = faceCentre(facei);
    }

    return tpts;
}


Foam::point Foam::boundBox::faceCentre(const direction facei) const
{
    point pt = boundBox::centre();

    switch (facei)
    {
        case 0: pt.x() = min().x(); break;  // 0: x-min, left
        case 1: pt.x() = max().x(); break;  // 1: x-max, right
        case 2: pt.y() = min().y(); break;  // 2: y-min, bottom
        case 3: pt.y() = max().y(); break;  // 3: y-max, top
        case 4: pt.z() = min().z(); break;  // 4: z-min, back
        case 5: pt.z() = max().z(); break;  // 5: z-max, front
        default:
        {
            FatalErrorInFunction
                << "Face:" << int(facei) << " should be [0..5]"
                << abort(FatalError);
        }
    }

    return pt;
}


void Foam::boundBox::reduce()
{
    Foam::reduce(min_, minOp<point>());
    Foam::reduce(max_, maxOp<point>());
}


Foam::boundBox Foam::boundBox::returnReduce(const boundBox& bb)
{
    boundBox work(bb);
    work.reduce();
    return work;
}


bool Foam::boundBox::intersects(const plane& pln) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    // Check as below(1) or above(2) - stop when it cuts both
    int side = 0;

    #undef  doLocalCode
    #define doLocalCode(Idx)                                                  \
    {                                                                         \
        side |= (pln.whichSide(hexCorner<Idx>()) == plane::BACK ? 1 : 2);     \
        if (side == 3) return true;  /* Both below and above: done */         \
    }

    // Each box corner
    doLocalCode(0);
    doLocalCode(1);
    doLocalCode(2);
    doLocalCode(3);
    doLocalCode(4);
    doLocalCode(5);
    doLocalCode(6);
    doLocalCode(7);

    #undef doLocalCode

    return false;
}


bool Foam::boundBox::intersects(const triPointRef& tri) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    // Simplest check - if any points are inside
    if (contains(tri.a()) || contains(tri.b()) || contains(tri.c()))
    {
        return true;
    }


    // Extent of box points projected onto axis
    const auto project_box = []
    (
        const boundBox& bb,
        const vector& axis,
        scalarMinMax& extent
    ) -> void
    {
        extent.reset(axis & bb.hexCorner<0>());
        extent.add(axis & bb.hexCorner<1>());
        extent.add(axis & bb.hexCorner<2>());
        extent.add(axis & bb.hexCorner<3>());
        extent.add(axis & bb.hexCorner<4>());
        extent.add(axis & bb.hexCorner<5>());
        extent.add(axis & bb.hexCorner<6>());
        extent.add(axis & bb.hexCorner<7>());
    };


    // Use separating axis theorem to determine if triangle and
    // (axis-aligned) bounding box intersect.

    scalarMinMax tri_extent(0);
    scalarMinMax box_extent(0);
    const boundBox& bb = *this;

    // 1.
    // Test separating axis defined by the box normals
    // (project triangle points)
    // - do first (largely corresponds to normal bound box rejection test)
    //
    // No intersection if extent of projected triangle points are outside
    // of the box range

    {
        // vector::X
        tri_extent.reset(tri.a().x());
        tri_extent.add(tri.b().x());
        tri_extent.add(tri.c().x());

        box_extent.reset(bb.min().x(), bb.max().x());

        if (!tri_extent.overlaps(box_extent))
        {
            return false;
        }

        // vector::Y
        tri_extent.reset(tri.a().y());
        tri_extent.add(tri.b().y());
        tri_extent.add(tri.c().y());

        box_extent.reset(bb.min().y(), bb.max().y());

        if (!tri_extent.overlaps(box_extent))
        {
            return false;
        }

        // vector::Z
        tri_extent.reset(tri.a().z());
        tri_extent.add(tri.b().z());
        tri_extent.add(tri.c().z());

        box_extent.reset(bb.min().z(), bb.max().z());

        if (!tri_extent.overlaps(box_extent))
        {
            return false;
        }
    }


    // 2.
    // Test separating axis defined by the triangle normal
    // (project box points)
    // - can use area or unit normal since any scaling is applied to both
    //   sides of the comparison.
    // - by definition all triangle points lie in the plane defined by
    //   the normal. It doesn't matter which of the points we use to define
    //   the triangle offset (extent) when projected onto the triangle normal

    vector axis = tri.areaNormal();

    tri_extent.reset(axis & tri.a());
    project_box(bb, axis, box_extent);

    if (!tri_extent.overlaps(box_extent))
    {
        return false;
    }


    // 3.
    // Test separating axes defined by the triangle edges, which are the
    // cross product of the edge vectors and the box face normals

    for (const vector& edgeVec : { tri.vecA(), tri.vecB(), tri.vecC() })
    {
        for (direction faceDir = 0; faceDir < vector::nComponents; ++faceDir)
        {
            axis = Zero;
            axis[faceDir] = 1;

            axis = (edgeVec ^ axis);

            // project tri
            tri_extent.reset(axis & tri.a());
            tri_extent.add(axis & tri.b());
            tri_extent.add(axis & tri.c());

            project_box(bb, axis, box_extent);

            if (!tri_extent.overlaps(box_extent))
            {
                return false;
            }
        }
    }

    return true;
}


bool Foam::boundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (!contains(p))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (contains(p))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBox::nearest(const point& p) const
{
    // Clip the point to the range of the bounding box
    return point
    (
        Foam::min(Foam::max(p.x(), min_.x()), max_.x()),
        Foam::min(Foam::max(p.y(), min_.y()), max_.y()),
        Foam::min(Foam::max(p.z(), min_.z()), max_.z())
    );
}


void Foam::boundBox::inflate(Random& rndGen, const scalar factor)
{
    vector newSpan(span());

    // Make 3D
    const scalar minSpan = factor * Foam::mag(newSpan);

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        newSpan[dir] = Foam::max(newSpan[dir], minSpan);
    }

    min_ -= cmptMultiply(factor*rndGen.sample01<vector>(), newSpan);
    max_ += cmptMultiply(factor*rndGen.sample01<vector>(), newSpan);
}


void Foam::boundBox::inflate
(
    Random& rndGen,
    const scalar factor,
    const scalar delta
)
{
    inflate(rndGen, factor);
    grow(delta);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::boundBox::operator&=(const boundBox& bb)
{
    min_ = ::Foam::max(min_, bb.min_);
    max_ = ::Foam::min(max_, bb.max_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBox)
        );
    }
    else
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstreamOption::BINARY)
    {
        Detail::readContiguous<boundBox>
        (
            is,
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }
    else
    {
        is >> bb.min_ >> bb.max_;
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
