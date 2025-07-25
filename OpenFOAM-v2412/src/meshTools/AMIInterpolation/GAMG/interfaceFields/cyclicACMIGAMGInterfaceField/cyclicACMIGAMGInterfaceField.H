/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::cyclicACMIGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface for Arbitrarily Coupled Mesh Interface
    (ACMI) fields.

SourceFiles
    cyclicACMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIGAMGInterfaceField_H
#define cyclicACMIGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "cyclicACMIGAMGInterface.H"
#include "cyclicACMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cyclicACMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public cyclicACMILduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const cyclicACMIGAMGInterface& cyclicACMIInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


        // Sending and receiving (distributed AMI)

            //- Current range of send requests (non-blocking)
            mutable labelRange sendRequests_;

            //- Current range of recv requests (non-blocking)
            mutable labelRange recvRequests_;

            //- Scalar send buffers
            mutable PtrList<List<solveScalar>> scalarSendBufs_;

            //- Scalar receive buffers
            mutable PtrList<List<solveScalar>> scalarRecvBufs_;


    // Private Member Functions

        //- No copy construct
        cyclicACMIGAMGInterfaceField
        (
            const cyclicACMIGAMGInterfaceField&
        ) = delete;

        //- No copy assignment
        void operator=(const cyclicACMIGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("cyclicACMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        cyclicACMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        cyclicACMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );

        //- Construct from GAMG interface and Istream
        cyclicACMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            Istream& is
        );

        //- Construct from GAMG interface and local and remote fields
        cyclicACMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& local,
            const UPtrList<lduInterfaceField>& other
        );

        //- Construct by assembling and return a clone.
        virtual autoPtr<GAMGInterfaceField> clone
        (
            const GAMGInterface& GAMGCp,
            const UPtrList<lduInterfaceField>& other
        ) const
        {
            return autoPtr<GAMGInterfaceField>
            (
                new cyclicACMIGAMGInterfaceField
                (
                    GAMGCp,
                    *this,      // local field
                    other       // other fields
                )
            );
        }


    //- Destructor
    virtual ~cyclicACMIGAMGInterfaceField() = default;


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicACMIInterface_.size();
            }


        // Interface matrix update

            //- Are all (receive) data available?
            virtual bool ready() const;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        //- Cyclic interface functions

            //- Does the interface field perform the transformation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicACMIInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicACMIInterface_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
