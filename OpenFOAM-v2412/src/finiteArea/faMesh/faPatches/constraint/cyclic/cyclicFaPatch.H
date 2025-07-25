/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::cyclicFaPatch

Description
    Cyclic-plane patch.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    cyclicFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicFaPatch_H
#define Foam_cyclicFaPatch_H

#include "coupledFaPatch.H"
#include "cyclicLduInterface.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFaPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicFaPatch
:
    public coupledFaPatch,
    public cyclicLduInterface
{
    // Private Member Functions

        void calcTransforms();

protected:

    // Protected Static Data

        //- Relative tolerance (for geometric matching). Is factor of
        //  maximum edge length per face.
        static const scalar matchTol_;


    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch geodesic distance between P and N
        void makeLPN(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;

        //- Make non-orthogonality correction vectors
        void makeCorrectionVectors(vectorField& cv) const { cv = Zero; }


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from dictionary
        cyclicFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm,
            const word& patchType
        );


    //- Destructor
    virtual ~cyclicFaPatch() = default;


    // Member Functions

        // Access

            //- Does this side own the patch ?
            //- Yes: it contains both sets of faces
            virtual bool owner() const
            {
                return true;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                NotImplemented;
                return index();
            }

            //- Return processor number
            virtual const cyclicLduInterface& neighbPatch() const
            {
                NotImplemented;
                return *this;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return coupledFaPatch::forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return coupledFaPatch::reverseT();
            }

            //- Initialise the calculation of the patch geometry
            virtual void initGeometry(PstreamBuffers&);

            //- Calculate the patch geometry
            virtual void calcGeometry(PstreamBuffers&);

            //- Initialise the patches for moving points
            virtual void initMovePoints(PstreamBuffers&, const pointField&);

            //- Correct patches after moving points
            virtual void movePoints(PstreamBuffers&, const pointField&);

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //- the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;


            //- Return the values of the given internal data adjacent to
            //- the interface as a field using edgeFace mapping
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData,
                const labelUList& edgeFaces
            ) const;


            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;


            //- Return neighbour field using edgeCells mapping
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData,
                const labelUList& edgeCells
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
