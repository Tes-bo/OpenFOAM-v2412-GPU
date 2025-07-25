/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

Class
    Foam::coupledFaPatch

Author
    Zeljko Tukovic and Hrvoje Jasak

Description
    coupledFaPatch is an abstract base class for patches that couple regions
    of the computational domain e.g. cyclic, arbitrary interfaces, sliding
    interfaces and processor-processor links.

SourceFiles
    coupledFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coupledFaPatch_H
#define Foam_coupledFaPatch_H

#include "lduInterface.H"
#include "faPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupledFaPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFaPatch
:
    public lduInterface,
    public faPatch
{
    // Private Data

        //- offset (distance) vector from one side of the couple to the other
        mutable vectorField separation_;

        //- Face transformation tensor
        mutable tensorField forwardT_;

        //- Neighbour-cell transformation tensor
        mutable tensorField reverseT_;


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const = 0;

        //- Make patch geodesic distance between P and N
        virtual void makeLPN(scalarField&) const = 0;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const = 0;

        //- Make non-orthogonality correction vectors
        virtual void makeCorrectionVectors(vectorField&) const = 0;

        //- Calculate the uniform transformation tensors
        void calcTransformTensors
        (
            const vector& Cf,
            const vector& Cr,
            const vector& nf,
            const vector& nr
        ) const;

        //- Calculate the transformation tensors
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr
        ) const;


public:

    //- Runtime type information
    TypeNameNoDebug("coupled");


    // Constructors

        //- Construct from components
        coupledFaPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label nbrPolyPatchIndex,
            const word& patchType
        )
        :
            faPatch(name, edgeLabels, index, bm, nbrPolyPatchIndex, patchType)
        {}

        //- Construct from dictionary
        coupledFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm,
            const word& patchType
        )
        :
            faPatch(name, dict, index, bm, patchType)
        {}


    //- Destructor
    virtual ~coupledFaPatch() = default;


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }

            //- Does this side own the patch ?
            virtual bool owner() const = 0;

            //- Does the coupled side own the patch ?
            virtual bool neighbour() const
            {
                return !owner();
            }

            //- The referring patchID (eg, for cyclics).
            //  Here for consistency with coupledPolyPatch
            virtual label referPatchID() const
            {
                return -1;
            }

            //- Are the coupled planes separated?
            bool separated() const
            {
                return separation_.size();
            }

            //- Return the offset (distance) vector from one side of the couple
            //  to the other
            const vectorField& separation() const
            {
                if (!separation_.size())
                {
                    FatalErrorInFunction
                        << "Coupled patches are not separated"
                        << abort(FatalError);
                }

                return separation_;
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                if (!forwardT_.size())
                {
                    FatalErrorInFunction
                        << "Coupled planes do not need transformation"
                        << abort(FatalError);
                }

                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                if (!reverseT_.size())
                {
                    FatalErrorInFunction
                        << "Coupled planes do not need transformation"
                        << abort(FatalError);
                }

                return reverseT_;
            }

            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return forwardT_.size() == 0;
            }


        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&) = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&) = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&) = 0;


        // Access functions for demand driven data

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const = 0;


        // Interface transfer functions

            //- Return faceCell addressing: lduInterface virtual function
            virtual const labelUList& faceCells() const
            {
                return edgeFaces();
            }

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const
            {}

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const = 0;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
