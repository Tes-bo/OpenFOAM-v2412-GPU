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
    Foam::pointMapper

Description
    This object provides mapping and fill-in information for point data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh.

SourceFiles
    pointMapper.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointMapper_H
#define Foam_pointMapper_H

#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                           Class pointMapper Declaration
\*---------------------------------------------------------------------------*/

class pointMapper
:
    public morphFieldMapper
{
    // Private Data

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- The size of the mapper = pointMesh::size()
        const label mapperLen_;

        //- Number of inserted (unmapped) points
        label nInsertedObjects_;

        //- Is the mapping direct
        bool direct_;


    // Demand-Driven Data

        //- Direct addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelList> directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable std::unique_ptr<labelListList> interpAddrPtr_;

        //- Interpolation weights
        mutable std::unique_ptr<scalarListList> weightsPtr_;

        //- Inserted points
        mutable std::unique_ptr<labelList> insertedObjectsPtr_;


    // Private Member Functions

        //- Calculate addressing for mapping with inserted points
        void calcAddressing() const;

public:

    // Generated Methods

        //- No copy construct
        pointMapper(const pointMapper&) = delete;

        //- No copy assignment
        void operator=(const pointMapper&) = delete;


    // Constructors

        //- Construct from mapPolyMesh
        pointMapper(const pointMesh&, const mapPolyMesh& mpm);



    //- Destructor
    virtual ~pointMapper();


    // Member Functions

        //- The mapper size
        virtual label size() const;

        //- Return size before mapping
        virtual label sizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        //- Are there unmapped values? i.e. do all size() elements get value
        virtual bool hasUnmapped() const
        {
            return insertedObjects();
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolation weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted points
        bool insertedObjects() const noexcept
        {
            return bool(nInsertedObjects_);
        }

        //- Return list of inserted points
        const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
