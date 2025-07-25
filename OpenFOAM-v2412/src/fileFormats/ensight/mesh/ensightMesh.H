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
    Foam::ensightMesh

Description
    Encapsulation of volume meshes for writing in ensight format.
    It manages cellZones, facesZone, patches.

    When cellZones are present (and not disabled), the cells are grouped
    in parts according to the zone.
    Any remaining \em unzoned cells are placed into the "internalMesh" part,
    which is always part 0. If cellZones are missing or disabled,
    all cells are placed into the "internalMesh" part.

    If one or more cellZones are explicitly requested, all other cells
    (including any unzoned cells) are ignored.

    The converted patch faces are restricted by the volume mesh coverage.
    Except when the entire internal mesh has been explicitly suppressed.

Note
    The internal data management uses a Map for cellZones, faceZones and
    patches. The internalMesh is treated as cellZone with a special index.

    Since the patches are subsetted by the internal mesh coverage,
    they are treated as indirect patches rather than regular poly patches.

SourceFiles
    ensightMesh.C
    ensightMeshI.H
    ensightMeshOptions.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightMesh_H
#define Foam_ensightMesh_H

#include "Map.H"
#include "ensightCells.H"
#include "ensightFaces.H"
#include "wordRes.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class ensightGeoFile;
class ensightMesh;

/*---------------------------------------------------------------------------*\
                         Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
public:

    // Forward Declarations
    class options;

    //- The zone-id for internal mesh or unzoned cells.
    static const label internalZone;


private:

    // Private Data

        //- Writer options
        const std::unique_ptr<options> options_;

        //- Reference to the OpenFOAM mesh
        const polyMesh& mesh_;

        //- Volume elements per cellZone, lookup by zone index.
        //  The zone -1 is reserved for internal mesh (unzoned cells)
        Map<ensightCells> cellZoneParts_;

        //- Face elements per faceZone, lookup by zone index.
        Map<ensightFaces> faceZoneParts_;

        //- Face elements per selected patch, lookup by patch index
        Map<ensightFaces> boundaryParts_;

        //- Track if it needs an update
        mutable bool needsUpdate_;

        //- Output verbosity level
        int verbose_;


    // Private Member Functions

        //- Clear all storage
        void clear();

        //- Enforce consistent index/part numbering
        void renumber();

        //- No copy construct
        ensightMesh(const ensightMesh&) = delete;

        //- No copy assignment
        void operator=(const ensightMesh&) = delete;


public:

    // Constructors

        //- Construct from mesh with all default options
        explicit ensightMesh(const polyMesh& mesh);

        //- Construct from components
        ensightMesh(const polyMesh& mesh, const options& opts);


    // Member Functions

        //- Output verbosity level
        int verbose() const noexcept;

        //- Change the output verbosity level.
        //  \return old level
        int verbose(const int level) noexcept;


    // Access

        //- Reference to the underlying polyMesh
        const polyMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Reference to the writer/mesh options
        inline const ensightMesh::options& option() const;

        //- Face elements per selected patch, lookup by patch index
        //  Process in sorted order.
        //  May require special treatment for zone -1 (internal).
        const Map<ensightCells>& cellZoneParts() const noexcept
        {
            return cellZoneParts_;
        }

        //- Face elements per faceZone, lookup by zone index.
        //  Process in sorted order.
        const Map<ensightFaces>& faceZoneParts() const noexcept
        {
            return faceZoneParts_;
        }

        //- Face elements per selected patch, lookup by patch index
        //  Process in sorted order.
        const Map<ensightFaces>& boundaryParts() const noexcept
        {
            return boundaryParts_;
        }


    // Sizing Information

        //- Any parts?
        inline bool empty() const noexcept;

        //- Number of parts
        inline label size() const noexcept;


    // Other

        //- Does the content need an update?
        bool needsUpdate() const noexcept
        {
            return needsUpdate_;
        }

        //- Mark as needing an update.
        //  May also free up unneeded data.
        //  Return false if already marked as expired.
        inline bool expire();

        //- Update for new mesh
        void correct();


    // Output

        //- Write geometry to file (normally in parallel).
        //- Adds beginGeometry() marker.
        //  If all geometry is disabled, it will simply writes the mesh
        //  bounding box (to ensure that the geometry file is non-empty)
        void write
        (
            ensightGeoFile& os,
            bool parallel = UPstream::parRun()
        ) const;

        //- Write geometry to file (normally in parallel).
        //- Adds beginGeometry() marker.
        FOAM_DEPRECATED_FOR(2024-05, "write(ensightGeoFile&, ...")
        inline void write
        (
            autoPtr<ensightGeoFile>& os,
            bool parallel = UPstream::parRun()
        ) const;
};


/*---------------------------------------------------------------------------*\
                    Class ensightMesh::options Declaration
\*---------------------------------------------------------------------------*/

//- Configuration options for the ensightMesh
class ensightMesh::options
{
    // Private Data

        //- Create in 'expired' mode
        bool lazy_;

        //- Use the internal mesh
        bool internal_;

        //- Use the boundary mesh
        bool boundary_;

        //- Handle cellZones (if internal_ is true)
        bool cellZones_;

        //- Selected patches only
        wordRes patchInclude_;

        //- Deselected patches
        wordRes patchExclude_;

        //- Selected cellZones
        wordRes cellZoneInclude_;

        //- Selected faceZones
        wordRes faceZoneInclude_;


public:

    // Constructors

        //- Default construct. Non-lazy with internal/boundary/cellZones.
        options();


    // Member Functions

    // Access

        //- Lazy creation? (ie, starts as needsUpdate)
        bool lazy() const noexcept;

        //- Using internal?
        bool useInternalMesh() const noexcept;

        //- Using boundary?
        bool useBoundaryMesh() const noexcept;

        //- Using faceZones?
        bool useFaceZones() const noexcept;

        //- Using cellZones?
        bool useCellZones() const noexcept;

        //- Selection of patches. Empty if unspecified.
        const wordRes& patchSelection() const noexcept
        {
            return patchInclude_;
        }

        //- Selection of black listed patches. Empty if unspecified.
        const wordRes& patchExclude() const noexcept
        {
            return patchExclude_;
        }

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& faceZoneSelection() const noexcept
        {
            return faceZoneInclude_;
        }

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& cellZoneSelection() const noexcept
        {
            return cellZoneInclude_;
        }


    // Edit

        //- Reset to defaults
        void reset();

        //- Lazy creation - ensightMesh starts as needsUpdate
        //  \return old value
        bool lazy(bool on) noexcept;

        //- Alter the useBoundaryMesh state
        //  \return old value
        bool useInternalMesh(bool on) noexcept;

        //- Alter the useBoundaryMesh state
        //  \return old value
        bool useBoundaryMesh(bool on);

        //- Alter the useCellZones state
        //  \return old value
        bool useCellZones(bool on);

        //- Define patch selection matcher
        void patchSelection(const UList<wordRe>& patterns);

        //- Define patch selection matcher
        void patchSelection(List<wordRe>&& patterns);

        //- Define patch selection to exclude
        void patchExclude(const UList<wordRe>& patterns);

        //- Define patch selection to exclude
        void patchExclude(List<wordRe>&& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(const UList<wordRe>& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(List<wordRe>&& patterns);

        //- Define cellZone selection matcher
        void cellZoneSelection(const UList<wordRe>& patterns);

        //- Define cellZone selection matcher
        void cellZoneSelection(List<wordRe>&& patterns);


    // Output

        //- Report values
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightMeshI.H"

#endif

// ************************************************************************* //
