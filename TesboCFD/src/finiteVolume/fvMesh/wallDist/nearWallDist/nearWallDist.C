/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020,2024 OpenCFD Ltd.
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

#include "nearWallDist.H"
#include "cellDistFuncs.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallDist::calculate()
{
    const cellDistFuncs wallUtils(mesh_);

    const volVectorField& cellCentres = mesh_.C();


    if (cellDistFuncs::useCombinedWallPatch)
    {
        // Collect indices of wall patches

        DynamicList<label> wallPatchIDs(mesh_.boundary().size());
        label nWalls = 0;

        forAll(mesh_.boundary(), patchi)
        {
            if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
            {
                wallPatchIDs.append(patchi);
                nWalls += mesh_.boundary()[patchi].size();
            }
            else
            {
                // Set distance to 0
                operator[](patchi) = 0.0;
            }
        }


        // Collect all mesh faces of wall patches

        DynamicList<label> faceLabels(nWalls);

        for (const label patchi : wallPatchIDs)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            const auto& pp = patch.patch();

            forAll(patch, i)
            {
                faceLabels.append(pp.start()+i);
            }
        }

        const uindirectPrimitivePatch wallPatch
        (
            UIndirectList<face>(mesh_.faces(), faceLabels),
            mesh_.points()
        );

        DynamicList<label> neighbours;

        nWalls = 0;
        for (const label patchi : wallPatchIDs)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            const labelUList& faceCells = patch.patch().faceCells();

            fvPatchScalarField& ypatch = operator[](patchi);

            forAll(patch, patchFacei)
            {
                // Get point connected neighbours (in wallPatch indices!)
                wallUtils.getPointNeighbours(wallPatch, nWalls, neighbours);

                label minFacei = -1;
                ypatch[patchFacei] = wallUtils.smallestDist
                (
                    cellCentres[faceCells[patchFacei]],
                    wallPatch,
                    neighbours,
                    minFacei
                );

                nWalls++;
            }
        }
    }
    else
    {
        // Get patch ids of walls
        const labelHashSet wallPatchIDs(wallUtils.getPatchIDs<wallPolyPatch>());

        // Size neighbours array for maximum possible
        DynamicList<label> neighbours(wallUtils.maxPatchSize(wallPatchIDs));


        // Correct all cells with face on wall

        forAll(mesh_.boundary(), patchi)
        {
            fvPatchScalarField& ypatch = operator[](patchi);

            const fvPatch& patch = mesh_.boundary()[patchi];

            if (isA<wallFvPatch>(patch))
            {
                const polyPatch& pPatch = patch.patch();

                const labelUList& faceCells = patch.faceCells();

                // Check cells with face on wall
                forAll(patch, patchFacei)
                {
                    wallUtils.getPointNeighbours
                    (
                        pPatch,
                        patchFacei,
                        neighbours
                    );

                    label minFacei = -1;

                    ypatch[patchFacei] = wallUtils.smallestDist
                    (
                        cellCentres[faceCells[patchFacei]],
                        pPatch,
                        neighbours,
                        minFacei
                    );
                }
            }
            else
            {
                ypatch = 0.0;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallDist::nearWallDist(const Foam::fvMesh& mesh)
:
    volScalarField::Boundary
    (
        mesh.boundary(),
        mesh.V(),           // Dummy internal field,
        fvPatchFieldBase::calculatedType()
    ),
    mesh_(mesh)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallDist::~nearWallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallDist::correct()
{
    if (mesh_.topoChanging())
    {
        const DimensionedField<scalar, volMesh>& V = mesh_.V();
        const fvBoundaryMesh& bnd = mesh_.boundary();

        this->setSize(bnd.size());
        forAll(*this, patchi)
        {
            this->set
            (
                patchi,
                fvPatchField<scalar>::New
                (
                    fvPatchFieldBase::calculatedType(),
                    bnd[patchi],
                    V
                )
            );
        }
    }

    calculate();
}


// ************************************************************************* //
