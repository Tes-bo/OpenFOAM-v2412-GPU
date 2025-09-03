/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::meshPointPatch

Description
    pointPatch with explicitly provided points instead of using the points
    of a polyPatch.

    Note: does not constrain displacement - is not a constraint patch.

SourceFiles
    meshPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef meshPointPatch_H
#define meshPointPatch_H

#include "pointPatch.H"
#include "polyPatch.H"
#include "autoPtr.H"
#include "patchIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class meshPointPatch Declaration
\*---------------------------------------------------------------------------*/

class meshPointPatch
:
    public pointPatch
{
private:

    // Private Member Functions

        //- No copy construct
        meshPointPatch(const meshPointPatch&) = delete;

        //- No copy assignment
        void operator=(const meshPointPatch&) = delete;


protected:

    // Protected Member Data

        //- Explicit mesh points
        const labelList meshPoints_;

        const List<pointConstraint> constraints_;

        //- Demand-driven local points
        mutable autoPtr<pointField> localPointsPtr_;

        //- Demand-driven local normals (assumes constructed with pointNormals
        //  or normal-only-constraint)
        mutable autoPtr<pointField> pointNormalsPtr_;


    // Protected Member Functions

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("meshPoint");


    // Constructors

        //- Construct from components
        meshPointPatch
        (
            const word& name,
            const labelUList& meshPoints,
            const List<pointConstraint>& constraints,
            const label index,
            const pointBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from single-constraint (i.e. slip, provided normals)
        meshPointPatch
        (
            const word& name,
            const labelUList& meshPoints,
            const vectorField& pointNormals,
            const label index,
            const pointBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        meshPointPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct given the original patch and a map
        meshPointPatch
        (
            const meshPointPatch& pp,
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
            return autoPtr<pointPatch>::NewFrom<meshPointPatch>
            (
                *this,
                bm,
                index,
                mapAddressing,
                reversePointMap
            );
        }


    //- Destructor
    virtual ~meshPointPatch() = default;


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return meshPoints().size();
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const
        {
            return meshPoints_;
        }

        //- Return constraints
        virtual const List<pointConstraint>& constraints() const
        {
            return constraints_;
        }

        //- Return pointField of points in patch
        virtual const pointField& localPoints() const;

        //- Return point unit normals. Assumes single constraint
        virtual const vectorField& pointNormals() const;

        //- Write the pointPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
