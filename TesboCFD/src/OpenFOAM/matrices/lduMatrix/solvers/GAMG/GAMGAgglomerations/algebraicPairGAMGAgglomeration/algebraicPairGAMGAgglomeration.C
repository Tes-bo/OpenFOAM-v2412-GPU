/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016,2023 OpenCFD Ltd.
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

#include "algebraicPairGAMGAgglomeration.H"
#include "lduMatrix.H"
#include "addToRunTimeSelectionTable.H"
//#include "cyclicAMILduInterface.H"
//#include "cyclicACMILduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(algebraicPairGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        algebraicPairGAMGAgglomeration,
        lduMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::algebraicPairGAMGAgglomeration::algebraicPairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict)
{
    if (matrix.hasLower())
    {
        agglomerate
        (
            nCellsInCoarsestLevel_,
            0,
            max(mag(matrix.upper()), mag(matrix.lower())),
            true
        );
    }
    else
    {
        agglomerate(nCellsInCoarsestLevel_, 0, mag(matrix.upper()), true);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::algebraicPairGAMGAgglomeration::movePoints()
{
    const bool ok = pairGAMGAgglomeration::movePoints();

    if (!requireUpdate_)
    {
        // TBD. For now commented out since cyclicAMI not known in this library.

/*
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
*/
    }

    return ok;
}


// ************************************************************************* //
