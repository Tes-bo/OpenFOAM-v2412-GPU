/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "lduPrimitiveMeshAssembly.H"

#include "cyclicFvPatch.H"
#include "cyclicAMIFvPatch.H"
#include "cyclicACMIFvPatch.H"

#include "lduPrimitiveProcessorInterface.H"
#include "AssemblyFvPatch.H"

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::lduPrimitiveMeshAssembly::update
(
    UPtrList<GeometricField<Type, fvPatchField, volMesh>>& psis
)
{
    label newFaces(0);
    label newFacesProc(0);
    label newPatches(0);

    const label nMeshes(meshes_.size());

    for (label i=0; i < nMeshes; ++i)
    {
        forAll(meshes_[i].interfaces(), patchI)
        {
            const polyPatch& pp = psis[i].mesh().boundaryMesh()[patchI];
            const fvPatchField<Type>& fvp = psis[i].boundaryField()[patchI];

            if (fvp.useImplicit())
            {
                if (pp.masterImplicit())
                {
                    label newFacesPatch(0);
                    label newFacesProcPatch(0);

                    pp.newInternalProcFaces(newFacesPatch, newFacesProcPatch);

                    newFaces += newFacesPatch;
                    newFacesProc += newFacesProcPatch;

                    if (newFacesProc > 0)
                    {
                        FatalErrorInFunction
                            << "The number of faces on either side of the "
                            << "coupled patch " << patchI << " are not "
                            << "the same. "
                            << "This might be due to the decomposition used. "
                            << "Please use decomposition preserving implicit "
                            << "patches on a single processor."
                            << exit(FatalError);
                    }

                    cellBoundMap_[i][patchI].setSize(newFacesPatch, -1);
                    facePatchFaceMap_[i][patchI].setSize(newFacesPatch, -1);
                    faceBoundMap_[i][patchI].setSize(newFacesPatch, -1);

                    const label nbrId = pp.neighbPolyPatchID();
                    const label meshNrbId = findNbrMeshId(pp, i);

                    cellBoundMap_[meshNrbId][nbrId].setSize
                    (
                        newFacesPatch,
                        -1
                    );
                    facePatchFaceMap_[meshNrbId][nbrId].setSize
                    (
                        newFacesPatch,
                        -1
                    );
                    faceBoundMap_[meshNrbId][nbrId].setSize
                    (
                        newFacesPatch,
                        -1
                    );
                }
            }
            else
            {
                patchMap_[i][patchI] = newPatches;
                patchLocalToGlobalMap_[i][patchI] = newPatches;
                newPatches++;
            }
        }
    }

    label virtualPatches = newPatches;

    // patchLocalToGlobalMap: map from original to asembled + extra Ids

    for (label i=0; i < nMeshes; ++i)
    {
        forAll(meshes_[i].interfaces(), patchI)
        {
            if (patchLocalToGlobalMap_[i][patchI] == -1)
            {
                patchLocalToGlobalMap_[i][patchI] = virtualPatches++;
            }
        }
    }

    DebugInfo << patchMap_ << endl;
    DebugInfo << patchLocalToGlobalMap_ << endl;

    label oldFaces(0);
    // Add the internal faces for each mesh
    for (label i=0; i < nMeshes; ++i)
    {
        newFaces += meshes_[i].lduAddr().upperAddr().size();
        oldFaces += meshes_[i].lduAddr().upperAddr().size();
    }

    if (debug)
    {
        Info<< " old total faces : " << oldFaces
            << " new total faces (internal) : " << newFaces
            << " new faces (remote) : " << newFacesProc
            << " new Faces : " << newFaces - oldFaces << endl;

        Info<< " new patches : " << newPatches << endl;

        Info<< "Local to Global patch Map  : "
            << patchLocalToGlobalMap_ << endl;
    }

    // This gives the global cellId given the local patchId for interfaces
    patchAddr_.setSize(newPatches);

    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes_[i].interfaces();

        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];

            if (globalPatchId != -1)
            {
                const labelUList& faceCells =
                    meshes_[i].lduAddr().patchAddr(patchI);

                // Fill local patchAddr for standard patches
                if (!faceCells.empty())
                {
                    patchAddr_[globalPatchId].setSize(faceCells.size(), -1);

                    for (label celli = 0; celli < faceCells.size(); ++celli)
                    {
                        patchAddr_[globalPatchId][celli] =
                            cellOffsets_[i] + faceCells[celli];
                    }
                }
            }
        }
    }

    // Interfaces
    interfaces().setSize(newPatches);
    // Primitive interfaces
    primitiveInterfaces().setSize(newPatches);

    // The interfaces are conserved (cyclics, proc, etc)
    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes_[i].interfaces();

        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];
            if (globalPatchId != -1)
            {
                if (interfacesLst.set(patchI))
                {
                    const polyPatch& pp =
                        psis[i].mesh().boundaryMesh()[patchI];

                    const fvBoundaryMesh& bm = psis[i].mesh().boundary();

                    if (isA<cyclicLduInterface>(interfacesLst[patchI]))
                    {
                        label nbrId = refCast
                            <const cyclicLduInterface>
                            (
                                interfacesLst[patchI]
                            ).neighbPatchID();

                        label globalNbr = patchMap()[i][nbrId];

                        primitiveInterfaces().set
                        (
                            globalPatchId,
                            new AssemblyFvPatch<cyclicFvPatch>
                            (
                                pp,
                                bm,
                                patchAddr_[globalNbr],
                                patchAddr_[globalPatchId],
                                globalNbr
                            )
                        );

                        interfaces().set
                        (
                            globalPatchId,
                            &primitiveInterfaces()[globalPatchId]
                        );
                    }
                    else if
                    (
                        isA<cyclicAMILduInterface>(interfacesLst[patchI])
                    )
                    {
                        label nbrId =
                            refCast
                            <
                                const cyclicAMIPolyPatch
                            >(pp).neighbPatchID();

                        label globalNbr = patchMap()[i][nbrId];

                        primitiveInterfaces().set
                        (
                            globalPatchId,
                            new AssemblyFvPatch<cyclicAMIFvPatch>
                            (
                                pp,
                                bm,
                                patchAddr_[globalNbr],
                                patchAddr_[globalPatchId],
                                globalNbr
                            )
                        );
                        interfaces().set
                        (
                            globalPatchId,
                            &primitiveInterfaces()[globalPatchId]
                        );
                    }
                    else if
                    (
                        isA<cyclicACMILduInterface>(interfacesLst[patchI])
                    )
                    {
                        label nbrId =
                            refCast
                            <
                                const cyclicACMIPolyPatch
                            >(pp).neighbPatchID();

                        label globalNbr = patchMap()[i][nbrId];

                        label nonOverlId =
                            refCast
                            <
                                const cyclicACMIPolyPatch
                            >(pp).nonOverlapPatchID();

                        label globalnonOverlId = patchMap()[i][nonOverlId];

                        primitiveInterfaces().set
                        (
                            globalPatchId,
                            new AssemblyFvPatch<cyclicACMIFvPatch>
                            (
                                pp,
                                bm,
                                patchAddr_[globalNbr],
                                patchAddr_[globalPatchId],
                                globalNbr,
                                globalnonOverlId
                            )
                        );
                        interfaces().set
                        (
                            globalPatchId,
                            &primitiveInterfaces()[globalPatchId]
                        );
                    }
                    else
                    {
                        primitiveInterfaces().set
                        (
                            globalPatchId,
                            nullptr
                        );
                        interfaces().set
                        (
                            globalPatchId,
                            interfacesLst.get(patchI)
                        );
                    }
                }
            }
        }
    }

    // Create new lower/upper addressing
    lowerAddr().setSize(newFaces, -1);
    upperAddr().setSize(newFaces, -1);

    label startIndex = 0;

    for (label i=0; i < nMeshes; ++i)
    {
        faceMap_[i].setSize(meshes_[i].lduAddr().lowerAddr().size(), -1);

        const label nFaces = meshes_[i].lduAddr().upperAddr().size();

        // Add individual addresses
        SubList<label>(lowerAddr(), nFaces, startIndex) =
            meshes_[i].lduAddr().lowerAddr();

        SubList<label>(upperAddr(), nFaces, startIndex) =
            meshes_[i].lduAddr().upperAddr();

        // Offset cellsIDs to global cell addressing
        label localFacei = 0;

        for (label facei=startIndex; facei < startIndex + nFaces; ++facei)
        {
            lowerAddr()[facei] += cellOffsets_[i];
            upperAddr()[facei] += cellOffsets_[i];

            faceMap_[i][localFacei++] = facei;
        }

        startIndex += nFaces;
    }
    // Add new lower/upper adressing for new internal faces corresponding
    // to patch faces that has a correspondent on the slave patch
    // (i.e map, AMI,etc)
    // Don't include faces that are in different proc

    label nFaces = startIndex;

    for (label i=0; i < nMeshes; ++i)
    {
        const lduInterfacePtrsList interfacesLst = meshes_[i].interfaces();

        forAll(interfacesLst, patchI)
        {
            const polyPatch& pp = psis[i].mesh().boundaryMesh()[patchI];

            const fvPatchField<Type>& fvp = psis[i].boundaryField()[patchI];

            if (fvp.useImplicit())
            {
                label meshNrbId = this->findNbrMeshId(pp, i);

                if (pp.masterImplicit())
                {
                    const labelUList& nbrFaceCells = pp.nbrCells();
                    const label nbrPatchId = pp.neighbPolyPatchID();
                    refPtr<labelListList> tlocalFaceToFace =
                        pp.mapCollocatedFaces();

                    const labelListList& localFaceToFace = tlocalFaceToFace();

                    // Compact target map
                    // key() = local face in proci
                    // *iter = compactId

                    label subFaceI(0);
                    forAll(pp.faceCells(), faceI)
                    {
                        const label cellI =
                            pp.faceCells()[faceI] + cellOffsets_[i];

                        const labelList& facesIds = localFaceToFace[faceI];

                        forAll(facesIds, j)
                        {
                            label nbrFaceId = facesIds[j];

                            // local faces
                            const label nbrCellI =
                                nbrFaceCells[nbrFaceId]
                              + cellOffsets_[meshNrbId];

                            lowerAddr()[nFaces] = min(cellI, nbrCellI);
                            upperAddr()[nFaces] = max(cellI, nbrCellI);

                            cellBoundMap_[i][patchI][subFaceI] = nbrCellI;
                            cellBoundMap_[meshNrbId][nbrPatchId][subFaceI] =
                                cellI;

                            facePatchFaceMap_[i][patchI][subFaceI] = faceI;
                            facePatchFaceMap_[meshNrbId][nbrPatchId][subFaceI]
                                = nbrFaceId;

                            faceBoundMap_[i][patchI][subFaceI] = nFaces;
                            faceBoundMap_[meshNrbId][nbrPatchId][subFaceI] =
                                nFaces;


                            ++subFaceI;
                            ++nFaces;
                        }
                    }
                }
            }
        }
    }

    if (newFaces != nFaces)
    {
       FatalErrorInFunction
            << "Incorrrect total number of faces in the assembled lduMatrix: "
            << newFaces << " != " << nFaces << nl
            << exit(FatalError);
    }

    // Sort upper-tri order
    {
        labelList oldToNew
        (
            upperTriOrder
            (
                lduAddr().size(),
                lowerAddr(),
                upperAddr()
            )
        );
        inplaceReorder(oldToNew, lowerAddr());
        inplaceReorder(oldToNew, upperAddr());

        for (labelList& faceMap : faceMap_)
        {
            for (label& facei : faceMap)
            {
                facei = oldToNew[facei];
            }
        }

        for (labelListList& bMap : faceBoundMap_)
        {
            for (labelList& faceMap : bMap)
            {
                for (label& facei : faceMap)
                {
                    if (facei != -1)
                    {
                        facei = oldToNew[facei];
                    }
                }
            }
        }
    }

    if (debug & 2)
    {
        DebugVar(faceBoundMap_);
        DebugVar(cellBoundMap_);
        DebugVar(lowerAddr());
        DebugVar(upperAddr());
        DebugVar(patchAddr_);
        DebugVar(cellOffsets_);
        DebugVar(faceMap_);
        (checkUpperTriangular(lduAddr().size(), lowerAddr(), upperAddr()));
        DebugVar(lduAddr().size())
    }
}


// ************************************************************************* //
