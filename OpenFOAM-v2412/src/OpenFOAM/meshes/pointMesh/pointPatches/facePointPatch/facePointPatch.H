/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::facePointPatch

Description
    A pointPatch based on a polyPatch

SourceFiles
    facePointPatch.C
    facePointPatchNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_facePointPatch_H
#define Foam_facePointPatch_H

#include "pointPatch.H"
#include "polyPatch.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class processorPointPatch;
class cyclicPointPatch;

/*---------------------------------------------------------------------------*\
                       Class facePointPatch Declaration
\*---------------------------------------------------------------------------*/

class facePointPatch
:
    public pointPatch
{
protected:

    // Protected Data

        //- Reference to the underlying polyPatch
        const polyPatch& polyPatch_;


    // Protected Member Functions

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
        facePointPatch(const facePointPatch&) = delete;

        //- No copy assignment
        void operator=(const facePointPatch&) = delete;


public:

    // Declare friendship with the coupledPointPatches to allow them to extend
    // the set of points with those not associated with faces
    friend class processorPointPatch;
    friend class cyclicPointPatch;


    //- Runtime type information
    TypeName(polyPatch::typeName_());


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            facePointPatch,
            polyPatch,
            (const polyPatch& patch, const pointBoundaryMesh& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch and boundary mesh
        facePointPatch
        (
            const polyPatch& p,
            const pointBoundaryMesh& pm
        );

        //- Construct given the original patch and a map
        facePointPatch
        (
            const facePointPatch& pp,
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
            return autoPtr<pointPatch>::NewFrom<facePointPatch>
            (
                *this,
                bm,
                index,
                mapAddressing,
                reversePointMap
            );
        }


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<facePointPatch> New
        (
            const polyPatch&,
            const pointBoundaryMesh&
        );


    //- Destructor
    virtual ~facePointPatch() = default;


    // Member Functions

        //- Return the polyPatch
        virtual const polyPatch& patch() const
        {
            return polyPatch_;
        }

        //- Return name
        virtual const word& name() const
        {
            return polyPatch_.name();
        }

        //- Return size
        virtual label size() const
        {
            return meshPoints().size();
        }

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const
        {
            return polyPatch_.index();
        }

        //- Return true if this patch field is coupled
        virtual bool coupled() const
        {
            return polyPatch_.coupled();
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const
        {
            return polyPatch_.meshPoints();
        }

        //- Return pointField of points in patch
        virtual const pointField& localPoints() const
        {
            return polyPatch_.localPoints();
        }

        //- Return point unit normals
        virtual const vectorField& pointNormals() const
        {
            return polyPatch_.pointNormals();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
