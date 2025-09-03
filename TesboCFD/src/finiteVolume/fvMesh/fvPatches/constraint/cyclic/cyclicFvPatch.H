/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::cyclicFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    cyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicFvPatch_H
#define Foam_cyclicFvPatch_H

#include "cyclicPolyPatch.H"
#include "coupledFvPatch.H"
#include "cyclicLduInterface.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicFvPatch
:
    public coupledFvPatch,
    public cyclicLduInterface
{
    // Private Data

        const cyclicPolyPatch& cyclicPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;


public:

    //- Runtime type information
    TypeName(cyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            cyclicPolyPatch_(refCast<const cyclicPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicPolyPatch& cyclicPatch() const
            {
                return cyclicPolyPatch_;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                return cyclicPolyPatch_.neighbPatchID();
            }

            virtual bool owner() const
            {
                return cyclicPolyPatch_.owner();
            }

            //- Return processor number
            virtual const cyclicFvPatch& neighbPatch() const
            {
                return refCast<const cyclicFvPatch>
                (
                    this->boundaryMesh()[cyclicPolyPatch_.neighbPatchID()]
                );
            }

            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return cyclicPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPolyPatch_.reverseT();
            }

            const cyclicFvPatch& neighbFvPatch() const
            {
                return refCast<const cyclicFvPatch>
                (
                    this->boundaryMesh()[cyclicPolyPatch_.neighbPatchID()]
                );
            }


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
            //- the interface as a field using a mapping faceCell
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData,
                const labelUList& faceCell
            ) const;


            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
