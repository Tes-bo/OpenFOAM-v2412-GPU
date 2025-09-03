/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::AssemblyFvPatch

Description
    An assembly of lduMatrix that is specific inter-region coupling
    through mapped patches.

SourceFiles
    AssemblyFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef AssemblyFvPatch_H
#define AssemblyFvPatch_H


#include "cyclicFvPatch.H"
#include "cyclicAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AssemblyFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class CyclicType>
class AssemblyFvPatch
:
    public CyclicType
{
protected:

    // Protected Data

        //- Nbr faceCell
        const labelList& nrbFaceCells_;

        //- Local faceCells
        const labelList& faceCells_;

        //- Nbr patch ID
        label nbrPatchID_;

        //- Non-overlap patch ID
        label nonOverlapPatchID_;


public:

    //- Runtime type information
    TypeName(CyclicType::typeName_());


    // Constructors

        //- Construct from components
        AssemblyFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm,
            const labelList& nrbFaceCells,
            const labelList& faceCells,
            const label nbrPatchID,
            const label nonOverPatchID = -1
        )
        :
            CyclicType(patch, bm),
            nrbFaceCells_(nrbFaceCells),
            faceCells_(faceCells),
            nbrPatchID_(nbrPatchID),
            nonOverlapPatchID_(nonOverPatchID)
        {}


    //- Destructor
    virtual ~AssemblyFvPatch() = default;


    // Member Functions

        //- Return neighbour field
        virtual tmp<labelField> internalFieldTransfer
        (
            const Pstream::commsTypes commsType,
            const labelUList& internalData
        ) const
        {
            auto tpif = tmp<labelField>::New(nrbFaceCells_.size());
            auto& pif = tpif.ref();

            forAll(pif, facei)
            {
                pif[facei] = internalData[nrbFaceCells_[facei]];
            }

            return tpif;
        }

        //- Return neighbour
        virtual label neighbPatchID() const
        {
            return nbrPatchID_;
        }

        //- Return non-overlap patch ID for ACMI
        virtual label nonOverlapPatchID() const
        {
            return nonOverlapPatchID_;
        }

        //- Return non-overlapping fvPatch for ACMI
        virtual const fvPatch& nonOverlapPatch() const
        {
            return this->boundaryMesh()[nonOverlapPatchID_];
        }

        //- Return nbr patch
        virtual const CyclicType& neighbPatch() const
        {
            return refCast<const CyclicType>
            (
                this->boundaryMesh()[nbrPatchID_]
            );
        }

        //- Return faceCells
        virtual const labelUList& faceCells() const
        {
            return faceCells_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeAssemblyFvPatch(Type)                                             \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(AssemblyFvPatch<Type>, 0);            \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
