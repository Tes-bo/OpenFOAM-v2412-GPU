/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::pointBoundaryMesh

Description
    Foam::pointBoundaryMesh

SourceFiles
    pointBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointBoundaryMesh_H
#define Foam_pointBoundaryMesh_H

#include "pointPatch.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class pointMesh;
class polyBoundaryMesh;
class wordRes;

/*---------------------------------------------------------------------------*\
                       Class pointBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class pointBoundaryMesh
:
    public pointPatchList,
    public regIOobject
{
    // Private Data

        //- Reference to mesh
        const pointMesh& mesh_;

        //- Demand-driven: list of patch ids per group
        mutable autoPtr<HashTable<labelList>> groupIDsPtr_;


    // Private Member Functions

        //- Calculate geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- Some patches have inGroup entries
        bool hasGroupIDs() const;

        //- Calculate group name to patch ids lookup
        void calcGroupIDs() const;

        //- Assign facePointPatches corresponding to the given polyBoundaryMesh
        void addPatches(const polyBoundaryMesh& pbm);

        //- No copy construct
        pointBoundaryMesh(const pointBoundaryMesh&) = delete;

        //- No copy assignment
        void operator=(const pointBoundaryMesh&) = delete;


public:

    //- Declare friendship with pointMesh
    friend class pointMesh;


    //- Runtime type information
    TypeName("pointBoundaryMesh");


    // Constructors

        //- Construct from polyBoundaryMesh
        pointBoundaryMesh(const pointMesh&, const polyBoundaryMesh&);

        //- Construct from IOobject and polyBoundaryMesh
        pointBoundaryMesh
        (
            const IOobject& io,
            const pointMesh&,
            const polyBoundaryMesh&
        );


    //- Destructor
    virtual ~pointBoundaryMesh() = default;


    // Member Functions

        //- Return the mesh reference
        const pointMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- The number of patches before the first processor patch.
        label nNonProcessor() const;

        //- The number of processorPointPatch patches
        label nProcessorPatches() const;

        //- Return a list of patch names
        wordList names() const;

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of physical types
        wordList physicalTypes() const;

        //- Return (sorted) patch indices for all matches.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRe& matcher, const bool useGroups) const;

        //- Return (sorted) patch indices for all matches.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRes& matcher, const bool useGroups) const;

        //- Return (sorted) patch indices for all selected matches that
        //- are not ignored.
        //- The selection logic as per Foam::wordRes::filter.
        //  Optionally matches patch groups.
        //  A no-op (returns empty list) for an empty select matcher
        labelList indices
        (
            const wordRes& select,
            const wordRes& ignore,
            const bool useGroups
        ) const;

        //- Find patch index given a name
        //  A no-op (returns -1) for an empty patchName
        label findPatchID
        (
            const word& patchName,
            const bool allowNotFound = true
        ) const;

        //- The patch indices per patch group
        const HashTable<labelList>& groupPatchIDs() const;

        //- Correct pointBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- Correct pointBoundaryMesh after topology update
        void updateMesh();

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle) If validBoundary calls updateMesh()
        //  after reordering to recalculate data (so call needs to be parallel
        //  sync in that case)
        void reorder(const labelUList& oldToNew, const bool validBoundary);

        //- writeData member function required by regIOobject
        virtual bool writeData(Ostream&) const;


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        FOAM_DEPRECATED_FOR(2018-08, "indices() method")
        labelList findIndices(const wordRe& key, bool useGroups) const
        {
            return indices(key, useGroups);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
