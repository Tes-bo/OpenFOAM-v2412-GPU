/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicAMIGAMGInterface.H"
#include "cyclicACMIGAMGInterface.H"
//#include "processorGAMGInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceAreaPairGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        faceAreaPairGAMGAgglomeration,
        lduMesh
    );

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        faceAreaPairGAMGAgglomeration,
        geometry
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceAreaPairGAMGAgglomeration::faceAreaPairGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(mesh, controlDict)
{
    const fvMesh& fvmesh = refCast<const fvMesh>(mesh);

    if (pairGAMGAgglomeration::requiresUpdate())
    {
        // This is the <=2406 logic. Apply perturbation to avoid jagged
        // agglomerations on axis-aligned meshes

        DebugPoutInFunction<< "Agglomerate with perturbation" << endl;

        //agglomerate(mesh, sqrt(fvmesh.magSf().primitiveField()));
        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,          //mesh,
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
    else
    {
        // In partial updating mode. Use Galilean invariant agglomeration
        // to e.g. have constant agglomeration for solid body rotation. Also
        // scale with face area to be consistent with (most) discretisation.

        DebugPoutInFunction<< "Agglomerate with faceArea magnitude" << endl;

        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,          //mesh,
            mag(fvmesh.magSf().primitiveField()),
            true
        );
    }
}


Foam::faceAreaPairGAMGAgglomeration::faceAreaPairGAMGAgglomeration
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(mesh, controlDict)
{
    if (pairGAMGAgglomeration::requiresUpdate())
    {
        // This is the <=2406 logic. Apply perturbation to avoid jagged
        // agglomerations on axis-aligned meshes

        DebugPoutInFunction<< "Agglomerate with perturbation" << endl;

        //agglomerate(mesh, sqrt(mag(faceAreas)));
        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,          //mesh,
            mag
            (
                cmptMultiply
                (
                    faceAreas
                   /sqrt(mag(faceAreas)),
                    vector(1, 1.01, 1.02)
                    //vector::one
                )
            ),
            true
        );
    }
    else
    {
        // In partial updating mode. Use Galilean invariant agglomeration
        // to e.g. have constant agglomeration for solid body rotation. Also
        // scale with face area to be consistent with (most) discretisation.

        DebugPoutInFunction<< "Agglomerate with faceArea magnitude" << endl;

        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,          //mesh,
            mag(faceAreas),
            true
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceAreaPairGAMGAgglomeration::movePoints()
{
    const bool ok = pairGAMGAgglomeration::movePoints();

    if (!requireUpdate_)
    {
        // movePoints lower down did not trigger update. Update in case of
        // cyclicAMI since contains addressing across patches and patches
        // have moved.

        bool hasCyclicAMI = false;
        if (!processorAgglomerate())
        {
            const auto& fineInterfaces = interfaceLevel(0);
            forAll(fineInterfaces, inti)
            {
                if (fineInterfaces.set(inti))
                {
                    const auto& intf = fineInterfaces[inti];
                    if
                    (
                        isA<cyclicAMILduInterface>(intf)
                     || isA<cyclicACMILduInterface>(intf)
                    )
                    {
                        hasCyclicAMI = true;

                        DebugPoutInFunction
                            << "Detected cyclicA(C)MI at interface " << inti
                            << ".Redoing patch agglomeration" << endl;

                        break;
                    }
                }
            }
        }


        if (hasCyclicAMI)
        {
            // Redo the interface agglomeration
            for
            (
                label fineLevelIndex = 0;
                fineLevelIndex < size();
                fineLevelIndex++
            )
            {
                // Near complete copy of boundary handling in
                // GAMGAgglomeration::agglomerateLduAddressing

                const auto& fineMesh = meshLevel(fineLevelIndex);
                const auto& fineInterfaces = interfaceLevel(fineLevelIndex);
                const lduAddressing& fineMeshAddr = fineMesh.lduAddr();

                // Get restriction map for current level
                const labelField& restrictMap =
                    restrictAddressing(fineLevelIndex);

                const label startOfRequests = UPstream::nRequests();
                forAll(fineInterfaces, inti)
                {
                    if (fineInterfaces.set(inti))
                    {
                        const auto& intf = fineInterfaces[inti];

                        if
                        (
                            isA<cyclicAMILduInterface>(intf)
                         || isA<cyclicACMILduInterface>(intf)
                        )
                        {
                            if (fineLevelIndex == 0)
                            {
                                intf.initInternalFieldTransfer
                                (
                                    Pstream::commsTypes::nonBlocking,
                                    restrictMap,
                                    fineMeshAddr.patchAddr(inti)
                                );
                            }
                            else
                            {
                                intf.initInternalFieldTransfer
                                (
                                    Pstream::commsTypes::nonBlocking,
                                    restrictMap
                                );
                            }
                        }
                    }
                }

                // Wait for comms to finised
                UPstream::waitRequests(startOfRequests);

                // New coarse-level interfaces
                //lduInterfacePtrsList coarseInterfaces(fineInterfaces.size());

                forAll(fineInterfaces, inti)
                {
                    if (fineInterfaces.set(inti))
                    {
                        const auto& intf = fineInterfaces[inti];

                        if
                        (
                            isA<cyclicAMILduInterface>(intf)
                         || isA<cyclicACMILduInterface>(intf)
                        )
                        {
                            tmp<labelField> restrictMapInternalField;

                            // The finest mesh uses patchAddr from the
                            // original lduAdressing.
                            // the coarser levels create thei own adressing for
                            // faceCells
                            if (fineLevelIndex == 0)
                            {
                                restrictMapInternalField =
                                    intf.interfaceInternalField
                                    (
                                        restrictMap,
                                        fineMeshAddr.patchAddr(inti)
                                    );
                            }
                            else
                            {
                                restrictMapInternalField =
                                    intf.interfaceInternalField
                                    (
                                        restrictMap
                                    );
                            }

                            tmp<labelField> nbrRestrictMapInternalField =
                                intf.internalFieldTransfer
                                (
                                    Pstream::commsTypes::nonBlocking,
                                    restrictMap
                                );

                            lduPrimitiveMesh& coarseMesh =
                                meshLevels_[fineLevelIndex];

                            autoPtr<GAMGInterface> coarseIntf
                            (
                                GAMGInterface::New
                                (
                                    inti,
                                    coarseMesh.rawInterfaces(),
                                    intf,
                                    restrictMapInternalField(),
                                    nbrRestrictMapInternalField(),
                                    fineLevelIndex,
                                    fineMesh.comm()
                                )
                            );

                            //coarseInterfaces.set(inti, coarseIntf.ptr());
                            coarseMesh.interfaces().set
                            (
                                inti,
                                coarseIntf.ptr()
                            );
                            coarseMesh.primitiveInterfaces().set
                            (
                                inti,
                               &coarseMesh.interfaces()[inti]
                            );
                        }
                    }
                }

                //meshLevels_[fineLevelIndex].addInterfaces
                //(
                //    coarseInterfaces,
                //    lduPrimitiveMesh::nonBlockingSchedule
                //    <
                //        processorGAMGInterface
                //    >
                //    (
                //        coarseInterfaces
                //    )
                //);
            }
        }

        if (debug)
        {
            printLevels();
        }
    }

    return ok;
}


// ************************************************************************* //
