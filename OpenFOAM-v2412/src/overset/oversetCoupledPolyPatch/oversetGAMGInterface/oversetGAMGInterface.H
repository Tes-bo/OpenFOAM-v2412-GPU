/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021,2023 OpenFOAM Foundation
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

Class
    Foam::oversetGAMGInterface

Description
    GAMG agglomerated cyclic AMI interface.

SourceFiles
    oversetGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oversetGAMGInterface_H
#define oversetGAMGInterface_H

#include "GAMGInterface.H"
#include "oversetLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class oversetGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class oversetGAMGInterface
:
    public GAMGInterface,
    virtual public oversetLduInterface
{
    // Private Member Functions

        //- No copy construct
        oversetGAMGInterface(const oversetGAMGInterface&) = delete;

        //- No copy assignment
        void operator=(const oversetGAMGInterface&) = delete;


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from fine level interface,
        //- local and neighbour restrict addressing
        oversetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from Istream
        oversetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );

        //- Construct from fine level interface,
        //- local and neighbour restrict addressing
        oversetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelList& interfaceMap,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddresssing,
            const labelUList& faceOffsets,
            const lduInterfacePtrsList& allInterfaces
        );

        //- Construct by assembling and returning a clone.
        virtual autoPtr<GAMGInterface> clone
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const labelList& interfaceMap,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddresssing,
            const labelUList& faceOffsets,
            const lduInterfacePtrsList& allInterfaces,
            const label coarseComm,
            const label myProcNo,
            const labelList& procAgglomMap
        ) const
        {
            return autoPtr<GAMGInterface>
            (
                new oversetGAMGInterface
                (
                    index,
                    coarseInterfaces,
                    *this,
                    interfaceMap,
                    faceCells,
                    faceRestrictAddresssing,
                    faceOffsets,
                    allInterfaces
                )
            );
        }


    //- Destructor
    virtual ~oversetGAMGInterface() = default;


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
