/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::fvPatchMapper

Description
    Mapping class for a fvPatchField.

SourceFiles
    fvPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvPatchMapper_H
#define Foam_fvPatchMapper_H

#include "fvPatchFieldMapper.H"
#include "faceMapper.H"
#include "fvPatch.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvPatch;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class fvPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class fvPatchMapper
:
    public fvPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to face mapper
        const faceMapper& faceMap_;

        //- Size before mapping
        const label sizeBeforeMapping_;


    // Demand-Driven Data

        mutable bool hasUnmapped_;

        //- Direct addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelList> directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelListList> interpAddrPtr_;

        //- Interpolation weights
        mutable std::unique_ptr<scalarListList> weightsPtr_;


    // Private Member Functions

        //- No copy construct
        fvPatchMapper(const fvPatchMapper&) = delete;

        //- No copy assignment
        void operator=(const fvPatchMapper&) = delete;


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;


public:

    // Constructors

        //- Construct from mappers
        fvPatchMapper
        (
            const fvPatch& patch,
            const faceMapper& faceMap
        );


    //- Destructor
    virtual ~fvPatchMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return patch_.size();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return faceMap_.direct();
        }

        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolation weights
        virtual const scalarListList& weights() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
