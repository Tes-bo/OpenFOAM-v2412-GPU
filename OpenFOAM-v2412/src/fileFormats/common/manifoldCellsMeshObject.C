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

#include "manifoldCellsMeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(manifoldCellsMeshObject);
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Change entries in labelList
static inline void replaceAll
(
    const label oldVal,
    const label newVal,
    labelUList& list
)
{
    for (label& val : list)
    {
        if (val == oldVal)
        {
            val = newVal;
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::refPtr<Foam::cellList> Foam::manifoldCellsMeshObject::filter
(
    const polyMesh& mesh,
    label& nCellsCorrected
)
{
    const auto& oldCells = mesh.cells();

    // Start with copy of existing cell list
    auto tnewCells = refPtr<cellList>::New(oldCells);
    auto& newCells = tnewCells.ref();

    DynamicList<label> removed;

    nCellsCorrected = 0;
    forAll(oldCells, celli)
    {
        const auto& oldCFaces = oldCells[celli];
        auto& newCFaces = newCells[celli];

        removed.clear();

        forAll(oldCFaces, cFacei)
        {
            const label facei = oldCFaces[cFacei];
            const label masteri = newCFaces[cFacei];

            const face& f = mesh.faces()[facei];

            forAll(oldCFaces, cFacej)
            {
                const label facej = oldCFaces[cFacej];
                const label masterj = newCFaces[cFacej];

                if (facej != facei)
                {
                    if (face::sameVertices(f, mesh.faces()[facej]))
                    {
                        //Pout<< "Same face:" << facei
                        //    << " master:" << masteri
                        //    << " verts:" << f << nl
                        //    << "         :" << facej
                        //    << " master:" << masterj
                        //    << " verts:" << mesh.faces()[facej]
                        //    << endl;

                        if (masteri < masterj)
                        {
                            replaceAll(masterj, masteri, newCFaces);
                            removed.append(masterj);
                        }
                        else if (masterj < masteri)
                        {
                            replaceAll(masteri, masterj, newCFaces);
                            removed.append(masteri);
                        }
                    }
                }
            }
        }

        if (removed.size())
        {
            // Compact out removed faces
            label newi = 0;
            for (const label facei : oldCFaces)
            {
                if (!removed.found(facei))
                {
                    newCFaces[newi++] = facei;
                }
            }
            newCFaces.resize(newi);
            ++nCellsCorrected;
        }
    }

    // Not needed (locally)
    if (nCellsCorrected == 0)
    {
        // Just use the existing mesh
        tnewCells.cref(mesh.cells());
    }

    // Number of cells corrected (globally)
    reduce(nCellsCorrected, sumOp<label>());

    return tnewCells;
}


Foam::refPtr<Foam::cellList> Foam::manifoldCellsMeshObject::filter
(
    const polyMesh& mesh
)
{
    label count = 0;
    return filter(mesh, count);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manifoldCellsMeshObject::manifoldCellsMeshObject(const polyMesh& mesh)
:
    MeshObject_type(mesh),
    cellsPtr_(nullptr),
    nCorrected_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::manifoldCellsMeshObject::manifold() const
{
    if (nCorrected_ < 0)
    {
        cellsPtr_ = filter(mesh(), nCorrected_);
    }

    return (nCorrected_ > 0);
}


const Foam::cellList& Foam::manifoldCellsMeshObject::cells() const
{
    if (nCorrected_ < 0)
    {
        cellsPtr_ = filter(mesh(), nCorrected_);
    }

    return (cellsPtr_ ? cellsPtr_.cref() : mesh().cells());
}


void Foam::manifoldCellsMeshObject::updateMesh(const mapPolyMesh&)
{
    cellsPtr_.reset(nullptr);
    nCorrected_ = -1;
}


// ************************************************************************* //
