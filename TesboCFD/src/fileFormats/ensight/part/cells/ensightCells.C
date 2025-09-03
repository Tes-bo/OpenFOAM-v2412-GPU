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

#include "ensightCells.H"
#include "error.H"
#include "bitSet.H"
#include "polyMesh.H"
#include "cellModel.H"
#include "manifoldCellsMeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightCells, 0);
}

const char* Foam::ensightCells::elemNames[5] =
{
    "tetra4", "pyramid5", "penta6", "hexa8", "nfaced"
};

static_assert
(
    Foam::ensightCells::nTypes == 5,
    "Support exactly 5 cell types (tetra4, pyramid5, penta6, hexa8, nfaced)"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightCells::resizeAll()
{
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

Foam::ensightCells::ensightCells()
:
    ensightPart(),
    manifold_(false),
    offsets_(Zero),
    sizes_(Zero)
{}


Foam::ensightCells::ensightCells(const string& description)
:
    ensightCells()
{
    rename(description);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::label, 5> Foam::ensightCells::sizes() const
{
    FixedList<label, 5> count;

    for (int typei = 0; typei < nTypes; ++typei)
    {
        count[typei] = size(elemType(typei));
    }

    return count;
}


Foam::label Foam::ensightCells::totalSize() const noexcept
{
    label count = 0;
    for (label n : sizes_)
    {
        count += n;
    }
    return count;
}


void Foam::ensightCells::clear()
{
    clearOut();

    ensightPart::clear();

    manifold_ = false;
    sizes_ = Zero;
    offsets_ = Zero;
}


void Foam::ensightCells::clearOut()
{}


void Foam::ensightCells::reduce()
{
    for (int typei = 0; typei < nTypes; ++typei)
    {
        sizes_[typei] = size(elemType(typei));
    }
    Foam::reduce(sizes_, sumOp<label>());
}


void Foam::ensightCells::sort()
{
    for (int typei=0; typei < nTypes; ++typei)
    {
        const labelRange sub(range(elemType(typei)));

        if (!sub.empty())
        {
            SubList<label> ids(addressing(), sub);

            Foam::sort(ids);
        }
    }
}


template<class Addressing>
void Foam::ensightCells::classifyImpl
(
    const polyMesh& mesh,
    const Addressing& cellIds
)
{
    manifold_ = manifoldCellsMeshObject::New(mesh).manifold();

    // References to cell shape models
    const cellModel& tet   = cellModel::ref(cellModel::TET);
    const cellModel& pyr   = cellModel::ref(cellModel::PYR);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& hex   = cellModel::ref(cellModel::HEX);

    const cellShapeList& shapes = mesh.cellShapes();

    // Pass 1: Count the shapes

    sizes_ = Zero;  // reset sizes
    for (const label id : cellIds)
    {
        const cellModel& model = shapes[id].model();

        elemType etype(elemType::NFACED);
        if (model == tet)
        {
            etype = elemType::TETRA4;
        }
        else if (model == pyr)
        {
            etype = elemType::PYRAMID5;
        }
        else if (model == prism)
        {
            etype = elemType::PENTA6;
        }
        else if (model == hex)
        {
            etype = elemType::HEXA8;
        }

        ++sizes_[etype];
    }

    resizeAll();    // adjust allocation
    sizes_ = Zero;  // reset sizes - use for local indexing here


    // Pass 2: Assign cell-id per shape type

    for (const label id : cellIds)
    {
        const cellModel& model = shapes[id].model();

        elemType etype(elemType::NFACED);
        if (model == tet)
        {
            etype = elemType::TETRA4;
        }
        else if (model == pyr)
        {
            etype = elemType::PYRAMID5;
        }
        else if (model == prism)
        {
            etype = elemType::PENTA6;
        }
        else if (model == hex)
        {
            etype = elemType::HEXA8;
        }

        add(etype, id);
    }
}


void Foam::ensightCells::classify(const polyMesh& mesh)
{
    // All mesh cells
    classifyImpl(mesh, labelRange(mesh.nCells()));
}


void Foam::ensightCells::classify
(
    const polyMesh& mesh,
    const labelUList& cellIds
)
{
    classifyImpl(mesh, cellIds);
}


void Foam::ensightCells::classify
(
    const polyMesh& mesh,
    const bitSet& selection
)
{
    classifyImpl(mesh, selection);
}


void Foam::ensightCells::writeDict(Ostream& os, const bool full) const
{
    os.beginBlock(type());

    os.writeEntry("id",     index()+1); // Ensight starts with 1
    os.writeEntry("name",   name());
    os.writeEntry("size",   size());

    if (full)
    {
        for (int typei=0; typei < ensightCells::nTypes; ++typei)
        {
            const auto etype = ensightCells::elemType(typei);

            os.writeKeyword(ensightCells::key(etype));

            cellIds(etype).writeList(os, 0) << endEntry;  // Flat output
        }
    }

    os.endBlock();
}


// ************************************************************************* //
