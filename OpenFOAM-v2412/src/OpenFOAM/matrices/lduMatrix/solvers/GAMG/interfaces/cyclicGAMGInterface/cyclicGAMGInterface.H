/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::cyclicGAMGInterface

Description
    GAMG agglomerated cyclic interface.

SourceFiles
    cyclicGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGAMGInterface_H
#define cyclicGAMGInterface_H

#include "GAMGInterface.H"
#include "cyclicLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicGAMGInterface
:
    public GAMGInterface,
    virtual public cyclicLduInterface
{
    // Private data

        //- Neighbour patch number
        label neighbPatchID_;

        //- Am I owner?
        bool owner_;

        //- Transformation tensor
        tensorField forwardT_;

        //- Transformation tensor
        tensorField reverseT_;


    // Private Member Functions

        //- No copy construct
        cyclicGAMGInterface(const cyclicGAMGInterface&) = delete;

        //- No copy assignment
        void operator=(const cyclicGAMGInterface&) = delete;


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        cyclicGAMGInterface
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
        cyclicGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );

        //- Construct from fine level interface,
        //- local and neighbour restrict addressing
        cyclicGAMGInterface
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
                new cyclicGAMGInterface
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
    virtual ~cyclicGAMGInterface() = default;


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Cyclic interface functions

            //- Return neighbour processor number
            virtual label neighbPatchID() const
            {
                return neighbPatchID_;
            }

            virtual bool owner() const
            {
                return owner_;
            }

            virtual const cyclicGAMGInterface& neighbPatch() const
            {
                return dynamic_cast<const cyclicGAMGInterface&>
                (
                    coarseInterfaces_[neighbPatchID()]
                );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return reverseT_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
