/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::ZoneMesh

Description
    A list of mesh zones.

SourceFiles
    ZoneMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ZoneMesh_H
#define Foam_ZoneMesh_H

#include "regIOobject.H"
#include "pointField.H"
#include "Map.H"
#include "HashSet.H"
#include "PtrList.H"
#include "bitSet.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class ZoneType, class MeshType> class ZoneMesh;

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream& os, const ZoneMesh<ZoneType, MeshType>& zones);


/*---------------------------------------------------------------------------*\
                           Class ZoneMesh Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType, class MeshType>
class ZoneMesh
:
    public PtrList<ZoneType>,
    public regIOobject
{
    // Private Data

        //- Reference to mesh
        const MeshType& mesh_;

        //- Demand-driven: map of zone labels for given element
        mutable std::unique_ptr<Map<label>> zoneMapPtr_;

        //- Demand-driven: map of additional zone labels for given element.
        //  Only populated if overlapping zones
        mutable std::unique_ptr<labelListList> additionalMapPtr_;

        //- Demand-driven: list of zone ids per group
        mutable std::unique_ptr<HashTable<labelList>> groupIDsPtr_;


    // Private Member Functions

        //- Total of number of addressed items (all zones)
        label totalSize() const;

        //- Create zone map
        void calcZoneMap() const;

        //- Some zones have inGroup entries
        bool hasGroupIDs() const;

        //- Calculate group name to zone ids lookup
        void calcGroupIDs() const;

        //- Clear addressing at this level
        void clearLocalAddressing();

        //- Populate/recreate from dictionary entries
        void populate(PtrList<entry>&& entries);

        //- Return true if contents were read
        //- (controlled by IOobject readOption flags).
        bool readIOcontents(const bool allowOptionalRead);


public:

    // Public Typedefs

        //- The zone type. Same as PtrList<ZoneType>::value_type
        typedef ZoneType zone_type;


    //- Debug switch to disallow the use of generic zones
    static int disallowGenericZones;


    // Generated Methods

        //- No copy construct
        ZoneMesh(const ZoneMesh&) = delete;

        //- No copy assignment
        void operator=(const ZoneMesh<ZoneType, MeshType>&) = delete;


    // Constructors

        //- Read construct from IOobject and mesh reference
        //- Any reading (mandatory, optional) based on IOobject properties.
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh
        );

        //- Construct empty with IOobject properties and a mesh reference.
        //- Does not read.
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            Foam::zero
        );

        //- Construct with specified size if not read.
        //- Any reading (mandatory, optional) based on IOobject properties.
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            const label size
        );

        //- Read construct (mandatory, optional) based on IOobject properties
        //- or use the fallback PtrList (with cloning).
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            const PtrList<ZoneType>& list
        );

        //- Read construct (mandatory, optional) based on IOobject properties
        //- or use the fallback PtrList (with cloning).
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            PtrList<entry>&& entries
        );


    //- Destructor
    ~ZoneMesh() = default;


    // Member Functions

        //- Return the mesh reference
        const MeshType& mesh() const noexcept { return mesh_; }

        //- Map of zones containing zone index for all zoned elements
        //  Return -1 if the object is not in the zone
        const Map<label>& zoneMap() const;

        //- Given a global object index, return the zone it is in.
        //  If object does not belong to any zones, return -1
        label whichZone(const label objectIndex) const;

        //- Given a global object index, return (in argument) its zones.
        //  Returns number of zones (0 if object does not belong to any zones)
        label whichZones
        (
            const label objectIndex,
            DynamicList<label>& zones
        ) const;

        //- Return a list of zone types
        wordList types() const;

        //- A list of the zone names
        wordList names() const;

        //- A list of the zone group names (if any)
        wordList groupNames() const;

        //- A list of zone names satisfying the input matcher
        wordList names(const wordRe& matcher) const;

        //- A list of zone names satisfying the input matchers
        wordList names(const wordRes& matcher) const;

        //- Sorted list of the zone names
        wordList sortedNames() const;

        //- Sorted list of zone names satisfying the input matcher
        wordList sortedNames(const wordRe& matcher) const;

        //- Sorted list of zone names satisfying the input matchers
        wordList sortedNames(const wordRes& matcher) const;


        //- Return (sorted) zone indices for all matches
        //  Optionally matches zone groups.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices
        (
            const wordRe& matcher,
            const bool useGroups = true
        ) const;

        //- Return (sorted) zone indices for all matches
        //  Optionally matches zone groups.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices
        (
            const wordRes& matcher,
            const bool useGroups = true
        ) const;

        //- Return (sorted) zone indices for all selected matches that
        //- are not ignored.
        //- The selection logic as per Foam::wordRes::filter.
        //  Optionally matches patch groups.
        //  A no-op (returns empty list) for an empty select matcher
        labelList indices
        (
            const wordRes& select,
            const wordRes& ignore,
            const bool useGroups = true
        ) const;

        //- Zone index for the first match, return -1 if not found
        //  A no-op (returns -1) for an empty key
        label findIndex(const wordRe& key) const;

        //- Zone index for the first match, return -1 if not found
        //  A no-op (returns -1) for an empty matcher
        label findIndex(const wordRes& matcher) const;

        //- Find zone index by name, return -1 if not found
        //  A no-op (returns -1) for an empty zoneName
        label findZoneID(const word& zoneName) const;

        //- Find zone by name and return const pointer, nullptr on error
        //  A no-op (returns nullptr) for an empty zoneName
        const ZoneType* cfindZone(const word& zoneName) const;

        //- Find zone by name and return pointer, nullptr on error
        //  A no-op (returns nullptr) for an empty zoneName
        ZoneType* findZone(const word& zoneName);


        //- Return all elements (cells, faces, points) contained in the
        //- listed zones.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        bitSet selection(const labelUList& zoneIds) const;

        //- Return all elements (cells, faces, points) that match the zone
        //- specification as a bitSet.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        //  Optionally matches zoneGroups.
        //  A no-op (returns empty bitSet) for an empty matcher
        bitSet selection
        (
            const wordRe& matcher,
            const bool useGroups = true
        ) const;

        //- Return all elements (cells, faces, points) that match the zone
        //- specification as a bitSet.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        //  A no-op (returns empty bitSet) for an empty matcher
        bitSet selection
        (
            const wordRes& matcher,
            const bool useGroups = true
        ) const;

        //- The zone indices per zone group
        const HashTable<labelList>& groupZoneIDs() const;

        //- Set/add group with zones
        void setGroup(const word& groupName, const labelUList& zoneIDs);

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have all zones and in same order.
        //  \return True if any errors.
        bool checkParallelSync(const bool report = false) const;

        //- Correct zone mesh after moving points
        void movePoints(const pointField& pts);


    // Storage Management

        //- Clear addressing
        void clearAddressing();

        //- Clear primitive addressing
        void clearPrimitives();

        //- Clear the zones
        void clear();

        //- The zoneMap has been allocated
        bool hasZoneMap() const noexcept { return bool(zoneMapPtr_); }


    // Member Operators

        //- Return const and non-const reference to zone by index.
        using PtrList<ZoneType>::operator[];

        //- Return const reference to zone by name.
        //  Fatal if the zone does not exist.
        const ZoneType& operator[](const word& zoneName) const;

        //- Return reference to an existing zone by name
        //  Fatal if the zone does not exist.
        ZoneType& operator[](const word& zoneName);

        //- Find an existing zone by name or create a new empty one
        //- if required.
        //
        //  To determine if the zone already existed or was newly created,
        //  it will be necessary to add additional logic in the caller.
        //  For example,
        //  \code
        //      const label nOrig = zones.size();
        //
        //      ZoneType& zn = zones("zoneName");
        //
        //      if (nOrig == zones.size()) { existing... } else { new... }
        //  \endcode
        //  \param zoneName the name of the zone
        //  \param verbose report if an existing zone was selected or
        //      a new zone was created.
        //  \return non-const reference to the existing or new zone
        ZoneType& operator()(const word& zoneName, const bool verbose=false);


    // IO

        //- Update internal meta-data (eg, prior to writing)
        void updateMetaData();

        //- The writeData member function required by regIOobject
        bool writeData(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <ZoneType, MeshType>
        (
            Ostream& os,
            const ZoneMesh<ZoneType, MeshType>& zones
        );


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        FOAM_DEPRECATED_FOR(2018-08, "indices() method")
        labelList findIndices(const wordRes& key) const
        {
            return indices(key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ZoneMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
