/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "ensightFaces.H"
#include "error.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightFaces, 0);
}

const char* Foam::ensightFaces::elemNames[3] =
{
    "tria3", "quad4", "nsided"
};

static_assert
(
    Foam::ensightFaces::nTypes == 3,
    "Support exactly 3 face types (tria3, quad4, nsided)"
);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Trivial shape classifier
inline Foam::ensightFaces::elemType whatType(const Foam::face& f) noexcept
{
    return
    (
        f.size() == 3
      ? Foam::ensightFaces::elemType::TRIA3
      : f.size() == 4
      ? Foam::ensightFaces::elemType::QUAD4
      : Foam::ensightFaces::elemType::NSIDED
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightFaces::resizeAll()
{
    // Invalidate any previous face ordering
    faceOrder_.clear();

    // Invalidate any previous flipMap
    flipMap_.clear();

    // Assign sub-list offsets, determine overall size

    label len = 0;

    auto iter = offsets_.begin();

    *iter = 0;
    for (const label n : sizes_)
    {
        len += n;

        *(++iter) = len;
    }

    // The addressing space
    addressing().resize(len, Zero);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFaces::ensightFaces()
:
    ensightPart(),
    faceOrder_(),
    flipMap_(),
    offsets_(Zero),
    sizes_(Zero)
{}


Foam::ensightFaces::ensightFaces(const string& description)
:
    ensightFaces()
{
    rename(description);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 3> Foam::ensightFaces::sizes() const
{
    FixedList<label, 3> count;

    for (int typei = 0; typei < nTypes; ++typei)
    {
        count[typei] = size(elemType(typei));
    }

    return count;
}


Foam::label Foam::ensightFaces::totalSize() const noexcept
{
    label count = 0;
    for (label n : sizes_)
    {
        count += n;
    }
    return count;
}


void Foam::ensightFaces::clear()
{
    clearOut();

    ensightPart::clear();

    faceOrder_.clear();
    flipMap_.clear();
    sizes_ = Zero;
    offsets_ = Zero;
}


void Foam::ensightFaces::clearOut()
{}


void Foam::ensightFaces::reduce()
{
    for (int typei = 0; typei < nTypes; ++typei)
    {
        sizes_[typei] = size(elemType(typei));
    }
    Foam::reduce(sizes_, sumOp<label>());
}


void Foam::ensightFaces::sort()
{
    // Some extra safety
    if (faceOrder_.size() != size())
    {
        faceOrder_.clear();
    }
    if (flipMap_.size() != size())
    {
        flipMap_.clear();
    }

    // Sort by face Ids.
    // Use to reorder flip maps and face-order too.

    for (int typei=0; typei < nTypes; ++typei)
    {
        const labelRange sub(range(elemType(typei)));

        if (!sub.empty())
        {
            SubList<label> ids(sub, addressing());
            labelList order(Foam::sortedOrder(ids));

            ids = reorder<labelList>(order, ids);

            // Sort flip map as well
            if (!flipMap_.empty())
            {
                SubList<bool> flips(flipMap_, sub);
                flips = reorder<boolList>(order, flips);
            }

            // Sort face ordering as well
            if (!faceOrder_.empty())
            {
                SubList<label> faceOrder(faceOrder_, sub);
                faceOrder = reorder<labelList>(order, faceOrder);
            }
        }
    }
}


void Foam::ensightFaces::classify(const UList<face>& faces)
{
    const label len = faces.size();

    // Pass 1: Count the shapes

    sizes_ = Zero;  // reset sizes
    for (label listi = 0; listi < len; ++listi)
    {
        const auto etype = whatType(faces[listi]);

        ++sizes_[etype];
    }

    resizeAll();    // adjust allocation
    sizes_ = Zero;  // reset sizes - use for local indexing here

    // Pass 2: Assign face-id per shape type

    for (label listi = 0; listi < len; ++listi)
    {
        const auto etype = whatType(faces[listi]);

        add(etype, listi);
    }
}


void Foam::ensightFaces::classify
(
    const UList<face>& faces,
    const labelRange& range
)
{
    const labelRange slice(range.subset0(faces.size()));

    // Operate on a local slice
    classify(SubList<face>(slice, faces));

    // Fixup to use the real faceIds instead of the 0-based slice
    incrAddressing(slice.start());
}


void Foam::ensightFaces::classify
(
    const UList<face>& faces,
    const labelUList& addr,
    const boolList& flipMap,
    const bitSet& exclude
)
{
    const label len = addr.size();
    const bool useFlip = (len == flipMap.size());

    // Pass 1: Count the shapes

    sizes_ = Zero;  // reset sizes
    for (label listi = 0; listi < len; ++listi)
    {
        const label faceId = addr[listi];

        if (!exclude.test(faceId))
        {
            const auto etype = whatType(faces[faceId]);

            ++sizes_[etype];
        }
    }

    resizeAll();    // adjust allocation
    sizes_ = Zero;  // reset sizes - use for local indexing here

    label nUsed = addressing().size();

    if (useFlip)
    {
        flipMap_.resize(nUsed);
        flipMap_ = false;
    }

    faceOrder_.resize(nUsed);

    // Pass 2: Assign face-id per shape type
    // - also record the face order

    nUsed = 0;
    for (label listi = 0; listi < len; ++listi)
    {
        const label faceId = addr[listi];

        if (!exclude.test(faceId))
        {
            const bool doFlip = useFlip && flipMap.test(listi);

            const auto etype = whatType(faces[faceId]);

            const label idx = add(etype, faceId, doFlip);

            faceOrder_[nUsed] = idx;
            ++nUsed;
        }
    }
}


void Foam::ensightFaces::writeDict(Ostream& os, const bool full) const
{
    os.beginBlock(type());

    os.writeEntry("id",     index()+1); // Ensight starts with 1
    os.writeEntry("name",   name());
    os.writeEntry("size",   size());

    if (full)
    {
        for (int typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const auto etype = ensightFaces::elemType(typei);

            os.writeKeyword(ensightFaces::key(etype));

            faceIds(etype).writeList(os, 0) << endEntry;  // Flat output
        }
    }

    os.endBlock();
}


// ************************************************************************* //
