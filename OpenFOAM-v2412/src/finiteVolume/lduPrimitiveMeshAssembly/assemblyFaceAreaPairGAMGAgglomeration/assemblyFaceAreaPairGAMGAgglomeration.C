/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "assemblyFaceAreaPairGAMGAgglomeration.H"
#include "lduMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "lduPrimitiveMeshAssembly.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(assemblyFaceAreaPairGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        assemblyFaceAreaPairGAMGAgglomeration,
        lduMatrix
    );

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        assemblyFaceAreaPairGAMGAgglomeration,
        geometry
    );
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

Foam::assemblyFaceAreaPairGAMGAgglomeration::
~assemblyFaceAreaPairGAMGAgglomeration()
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::assemblyFaceAreaPairGAMGAgglomeration::
assemblyFaceAreaPairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict)
{
    const lduMesh& ldumesh = matrix.mesh();

    if (isA<lduPrimitiveMeshAssembly>(ldumesh))
    {
        const lduPrimitiveMeshAssembly& mesh =
            refCast<const lduPrimitiveMeshAssembly>(ldumesh);

        vectorField faceAreas(mesh.lduAddr().upperAddr().size(), Zero);

        const labelListList& faceMap = mesh.faceMap();

        for (label i=0; i < mesh.meshes().size(); ++i)
        {
            const fvMesh& m = refCast<const fvMesh>(mesh.meshes()[i]);
            const labelList& subFaceMap = faceMap[i];
            const vectorField& areas = m.Sf().internalField();

            forAll(subFaceMap, facei)
            {
                faceAreas[subFaceMap[facei]] = areas[facei];
            }

            const polyBoundaryMesh& patches = m.boundaryMesh();

            // Fill faceAreas for new faces
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];
                label globalPatchID = mesh.patchMap()[i][patchI];

                if (globalPatchID == -1)
                {
                    if (pp.masterImplicit())
                    {
                        const vectorField& sf = m.boundary()[patchI].Sf();

                        if (isA<cyclicAMIPolyPatch>(pp))
                        {
                            const cyclicAMIPolyPatch& mpp =
                                    refCast<const cyclicAMIPolyPatch>(pp);

                            const scalarListList& srcWeight =
                                mpp.AMI().srcWeights();

                            label subFaceI = 0;
                            forAll(pp.faceCells(), faceI)
                            {
                                const scalarList& w = srcWeight[faceI];

                                for(label j=0; j<w.size(); j++)
                                {
                                    const label globalFaceI =
                                        mesh.faceBoundMap()[i][patchI][subFaceI];

                                    if (globalFaceI != -1)
                                    {
                                        faceAreas[globalFaceI] = w[j]*sf[faceI];
                                    }
                                    subFaceI++;
                                }
                            }
                        }
                        else if (isA<cyclicACMIPolyPatch>(pp))
                        {
                            const cyclicACMIPolyPatch& mpp =
                                refCast<const cyclicACMIPolyPatch>(pp);

                            const scalarListList& srcWeight =
                                mpp.AMI().srcWeights();

                            const scalarField& mask = mpp.mask();
                            const scalar tol = mpp.tolerance();
                            label subFaceI = 0;
                            forAll(mask, faceI)
                            {
                                const scalarList& w = srcWeight[faceI];

                                for(label j=0; j<w.size(); j++)
                                {
                                    if (mask[faceI] > tol)
                                    {
                                        const label globalFaceI =
                                            mesh.faceBoundMap()[i]
                                            [patchI][subFaceI];

                                        faceAreas[globalFaceI] = w[j]*sf[faceI];
                                    }
                                    subFaceI++;
                                }
                            }
                        }
                        else
                        {
                            forAll(pp.faceCells(), faceI)
                            {
                                const label globalFaceI =
                                    mesh.faceBoundMap()[i][patchI][faceI];

                                if (globalFaceI != -1)
                                {
                                    faceAreas[globalFaceI] = sf[faceI];
                                }
                            }
                        }
                    }
                }
            }
        }

        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,      //mesh,
            mag
            (
                cmptMultiply
                (
                    faceAreas/sqrt(mag(faceAreas)),
                    vector(1, 1.01, 1.02)
                )
            ),
            true
        );
    }
    else
    {
        // Revert to faceAreaPairGAMGAgglomeration
        const fvMesh& fvmesh = refCast<const fvMesh>(matrix.mesh());

        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,      //matrix.mesh(),
            mag
            (
                cmptMultiply
                (
                    fvmesh.Sf().primitiveField()
                   /sqrt(fvmesh.magSf().primitiveField()),
                    vector(1, 1.01, 1.02)
                    //vector::one
                )
            ),
            true
        );
    }
}


Foam::assemblyFaceAreaPairGAMGAgglomeration::
assemblyFaceAreaPairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict)
{
    agglomerate
    (
        nCellsInCoarsestLevel_,
        0,      //matrix.mesh(),
        mag
        (
            cmptMultiply
            (
                faceAreas/sqrt(mag(faceAreas)),
                vector(1, 1.01, 1.02)
            )
        ),
        true
    );
}


// ************************************************************************* //
