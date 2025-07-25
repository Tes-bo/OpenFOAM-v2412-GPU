/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::faAreaMapper

Description
    FA area mapper.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faAreaMapper.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faAreaMapper_H
#define Foam_faAreaMapper_H

#include "morphFieldMapper.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faAreaMapper Declaration
\*---------------------------------------------------------------------------*/

class faAreaMapper
:
    public morphFieldMapper
{
    // Private Data

        //- Reference to mesh mapper
        const faMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Size before mapping (old mesh size)
        const label sizeBeforeMapping_;

        //- Is the mapping direct
        bool direct_;


    // Demand-Driven Data

        mutable bool hasUnmapped_;

        //- New face labels after mapping
        mutable std::unique_ptr<labelList> newFaceLabelsPtr_;

        //- New face labels after mapping
        mutable std::unique_ptr<labelList> newFaceLabelsMapPtr_;


        //- Direct addressing (only one form of addressing is used)
        mutable std::unique_ptr<labelList> directAddrPtr_;

        //- Interpolated addressing (only one form of addressing is used)
        mutable std::unique_ptr<labelListList> interpAddrPtr_;

        //- Interpolation weights
        mutable std::unique_ptr<scalarListList> weightsPtr_;

        //- Inserted faces
        mutable std::unique_ptr<labelList> insertedObjectsPtr_;


    // Private Member Functions

        //- Calculate addressing
        void calcAddressing() const;


public:

    // Generated Methods

        //- No copy construct
        faAreaMapper(const faAreaMapper&) = delete;

        //- No copy assignment
        void operator=(const faAreaMapper&) = delete;


    // Constructors

        //- Construct from components
        faAreaMapper
        (
            const faMesh& mesh,
            const mapPolyMesh& mpm
        );


    //- Destructor
    virtual ~faAreaMapper();


    // Member Functions

        //- Return new face labels
        const labelList& newFaceLabels() const;

        //- Return new face labels map
        //  For new faces return old face index if it exists
        //  If the face has been added, index will be -1
        const labelList& newFaceLabelsMap() const;

        //- Return size
        virtual label size() const
        {
            return newFaceLabels().size();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
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

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return !insertedObjectLabels().empty();
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
