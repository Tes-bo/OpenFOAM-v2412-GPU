/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "oversetGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        oversetGAMGInterface,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        oversetGAMGInterface,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetGAMGInterface::oversetGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces
    )
{
    // Problem:
    // - we only want the oversetGAMGInterfaceField
    // - but this needs a oversetGAMGInterface of at least size 1
    //   (since faceRestrictAddressing cannot have -1 or so in it - it
    //    can only have >= 0 i.e. we cannot agglomerate into nothing)
    // - note that the end result is not used - the (current) corresponding
    //   oversetGAMGInterfaceField has dummy updateInterfaceMatrix.

    // Construct face agglomeration from cell agglomeration. Code same as
    // in e.g. cyclicAMIGAMGInterface.
    {
        // From coarse face to cell
        DynamicList<label> dynFaceCells(localRestrictAddressing.size());

        // From face to coarse face
        DynamicList<label> dynFaceRestrictAddressing
        (
            localRestrictAddressing.size()
        );

        Map<label> masterToCoarseFace(localRestrictAddressing.size());

        for (const label curMaster : localRestrictAddressing)
        {
            const auto iter = masterToCoarseFace.cfind(curMaster);

            if (iter.good())
            {
                // Already have coarse face
                dynFaceRestrictAddressing.append(iter.val());
            }
            else
            {
                // New coarse face
                const label coarseI = dynFaceCells.size();

                dynFaceRestrictAddressing.append(coarseI);
                dynFaceCells.append(curMaster);
                masterToCoarseFace.insert(curMaster, coarseI);
            }
        }

        faceCells_.transfer(dynFaceCells);
        faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
    }
}


Foam::oversetGAMGInterface::oversetGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    GAMGInterface(index, coarseInterfaces, is)
{}


Foam::oversetGAMGInterface::oversetGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelList& interfaceMap,
    const labelUList& faceCells,
    const labelUList& faceRestrictAddresssing,
    const labelUList& faceOffsets,
    const lduInterfacePtrsList& allInterfaces
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces,
        faceCells,
        faceRestrictAddresssing
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::oversetGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return tmp<labelField>::New(iF);
}


void Foam::oversetGAMGInterface::write(Ostream& os) const
{
    GAMGInterface::write(os);
}


// ************************************************************************* //
