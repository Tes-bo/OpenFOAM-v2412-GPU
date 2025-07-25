/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::oversetFvPatch

Description
    Patch for indicating interpolated boundaries (in overset meshes).

SourceFiles
    oversetFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatch_H
#define oversetFvPatch_H

#include "coupledFvPatch.H"
#include "oversetLduInterface.H"
#include "lduInterface.H"
#include "fvPatch.H"
#include "oversetPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oversetFvPatch Declaration
\*---------------------------------------------------------------------------*/

class oversetFvPatch
:
    public fvPatch,
    public lduInterface,
    public oversetLduInterface
{
    // Private Data

        const oversetPolyPatch& oversetPolyPatch_;


public:

    //- Runtime type information
    TypeName(oversetPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        oversetFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            oversetPolyPatch_(refCast<const oversetPolyPatch>(patch))
        {}


    // Member Functions

        // Access

            //- Return local reference cast into the overset patch
            const oversetPolyPatch& oversetPatch() const
            {
                return oversetPolyPatch_;
            }

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return oversetPolyPatch_.faceCells();
            }

            //- Am I the master interface
            virtual bool master() const
            {
                return oversetPolyPatch_.master();
            }


        // Interface Transfer Functions

            //- Return the values of the given internal data adjacent to
            //- the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Return the values of the given internal data adjacent to
            //- the interface as a field using faceCell mapping
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData,
                const labelUList& faceCells
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            //- using faceCells mapping
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF,
                const labelUList& faceCells
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            //- using nbrFaceCells mapping
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF,
                const labelUList& nbrFaceCells
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
