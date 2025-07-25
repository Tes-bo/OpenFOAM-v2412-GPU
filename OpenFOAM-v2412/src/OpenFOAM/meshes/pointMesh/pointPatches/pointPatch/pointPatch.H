/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::pointPatch

Description
    Basic pointPatch represents a set of points from the mesh.

SourceFiles
    pointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointPatch_H
#define Foam_pointPatch_H

#include "patchIdentifier.H"
#include "labelList.H"
#include "vectorField.H"
#include "pointField.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class pointBoundaryMesh;
class pointConstraint;
class pointPatch;
class PstreamBuffers;

Ostream& operator<<(Ostream&, const pointPatch&);

//- Store lists of pointPatch as a PtrList
typedef PtrList<pointPatch> pointPatchList;

/*---------------------------------------------------------------------------*\
                      Class pointPatch Declaration
\*---------------------------------------------------------------------------*/

class pointPatch
:
    public patchIdentifier
{
    // Private Data

        //- Reference to boundary mesh
        const pointBoundaryMesh& boundaryMesh_;


protected:

    // Protected Member Functions

        // The geometry initialisation is called by pointBoundaryMesh
        friend class pointBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&)
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&)
        {}

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&)
        {}


        //- No copy construct
        pointPatch(const pointPatch&) = delete;

        //- No copy assignment
        void operator=(const pointPatch&) = delete;


public:

    //- Runtime type information
    TypeName("basePatch");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatch,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const pointBoundaryMesh& bm,
                const word& patchType
            ),
            (name, dict, index, bm, patchType)
        );


    // Constructor

        //- Construct from components
        explicit pointPatch
        (
            const word& name,
            const label index,
            const pointBoundaryMesh& bm,
            const word& physicalType,
            const wordList& inGroups
        )
        :
            patchIdentifier(name, index, physicalType, inGroups),
            boundaryMesh_(bm)
        {}

        //- Construct from dictionary
        explicit pointPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointBoundaryMesh& bm
        )
        :
            patchIdentifier(name, dict, index),
            boundaryMesh_(bm)
        {}

        //- Construct given the original patch and a map
        explicit pointPatch
        (
            const pointPatch& pp,
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        )
        :
            patchIdentifier(pp.name(), index, pp.physicalType(), pp.inGroups()),
            boundaryMesh_(bm)
        {}

        //- Construct and return a subset clone,
        //- resetting the point list and boundary mesh
        virtual autoPtr<pointPatch> clone
        (
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        ) const = 0;


    // Selectors

        //- Return a pointer to a new patch created on freestore. Returns
        //- null if not found
        static autoPtr<pointPatch> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const pointBoundaryMesh&
        );


    //- Destructor
    virtual ~pointPatch() = default;


    // Member Functions

        //- Return size
        virtual label size() const = 0;

        //- Return boundaryMesh reference
        const pointBoundaryMesh& boundaryMesh() const
        {
            return boundaryMesh_;
        }

        //- Return true if this patch field is coupled
        virtual bool coupled() const
        {
            return false;
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const = 0;

        //- Return pointField of points in patch
        virtual const vectorField& localPoints() const = 0;

        //- Return point unit normals
        virtual const vectorField& pointNormals() const = 0;

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return word::null;
        }

        //- Accumulate the effect of constraint direction of this patch
        virtual void applyConstraint
        (
            const label pointi,
            pointConstraint&
        ) const
        {}

        //- Write the pointPatch data as a dictionary
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const pointPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
