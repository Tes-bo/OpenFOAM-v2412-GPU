/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::processorPointPatch

Description
    Processor patch boundary needs to be such that the ordering of
    points in the patch is the same on both sides.

    Looking at the creation of the faces on both sides of the processor
    patch they need to be identical on both sides with the normals pointing
    in opposite directions.  This is achieved by calling the reverseFace
    function in the decomposition.  It is therefore possible to re-create
    the ordering of patch points on the neighbour side by reversing all the
    patch faces of the owner.

SourceFiles
    processorPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorPointPatch_H
#define Foam_processorPointPatch_H

#include "coupledFacePointPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorPointPatch Declaration
\*---------------------------------------------------------------------------*/

class processorPointPatch
:
    public coupledFacePointPatch
{
    // Private Data

        const processorPolyPatch& procPolyPatch_;

        mutable labelList reverseMeshPoints_;


    // Private Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


        //- No copy construct
        processorPointPatch(const processorPointPatch&) = delete;

        //- No copy assignment
        void operator=(const processorPointPatch&) = delete;

public:

    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );

        //- Construct given the original patch and a map
        processorPointPatch
        (
            const processorPointPatch& patch,
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        );

        //- Construct and return a subset clone,
        //- resetting the point list and boundary mesh
        virtual autoPtr<pointPatch> clone
        (
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        ) const
        {
            return autoPtr<pointPatch>::NewFrom<processorPointPatch>
            (
                *this,
                bm,
                index,
                mapAddressing,
                reversePointMap
            );
        }


    //- Destructor
    virtual ~processorPointPatch() = default;


    // Member Functions

        //- Return message tag to use for communication
        virtual int tag() const
        {
            return procPolyPatch_.tag();
        }

        //- Return communicator used for communication
        virtual label comm() const
        {
            return procPolyPatch_.comm();
        }

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Return processor number
        int myProcNo() const
        {
            return procPolyPatch_.myProcNo();
        }

        //- Return neighbour processor number
        int neighbProcNo() const
        {
            return procPolyPatch_.neighbProcNo();
        }

        //- Does the processor own the patch ?
        bool owner() const
        {
            return procPolyPatch_.owner();
        }

        //- Is the processor the patch neighbour ?
        bool neighbour() const
        {
            return !owner();
        }

        //- Return the underlying processorPolyPatch
        const processorPolyPatch& procPolyPatch() const
        {
            return procPolyPatch_;
        }

        //- Return mesh points in the correct order for the receiving side
        const labelList& reverseMeshPoints() const
        {
            return reverseMeshPoints_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
