/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::triSurface

Description
    Triangulated surface description with patch information.

SourceFiles
    triSurface.C
    triSurfaceAddressing.C
    triSurfaceIO.C
    triSurfaceStitch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_triSurface_H
#define Foam_triSurface_H

#include "primitivePatch.H"
#include "PatchTools.H"
#include "labelledTri.H"
#include "boolList.H"
#include "bitSet.H"
#include "HashSet.H"
#include "geometricSurfacePatchList.H"
#include "surfacePatchList.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class IOobject;
class IFstream;
class surfZone;
class triSurface;

template<class Face> class MeshedSurface;

Istream& operator>>(Istream&, triSurface&);
Ostream& operator<<(Ostream&, const triSurface&);


/*---------------------------------------------------------------------------*\
                         Class triSurface Declaration
\*---------------------------------------------------------------------------*/

class triSurface
:
    public PrimitivePatch<::Foam::List<labelledTri>, pointField>
{
    // Private Typedefs

        //- Internal mesh storage type
        typedef PrimitivePatch<::Foam::List<labelledTri>, pointField>
            MeshReference;


    // Private Data

        //- Patch information
        //  (face ordering nFaces/startFace only used during reading, writing)
        geometricSurfacePatchList patches_;


    // Demand Driven

        //- Edge-face addressing (sorted)
        mutable std::unique_ptr<labelListList> sortedEdgeFacesPtr_;

        //- Label of face that 'owns' edge
        //- i.e. e.vec() is righthanded walk along face
        mutable std::unique_ptr<labelList> edgeOwnerPtr_;


    // Private Member Functions

        //- Calculate sorted edgeFaces
        void calcSortedEdgeFaces() const;

        //- Calculate owner
        void calcEdgeOwner() const;

        //- Sort faces according to region.
        //  Returns patch list and sets faceMap to index of labelledTri
        //  inside *this.
        surfacePatchList calcPatches(labelList& faceMap) const;

        //- Sets default values for patches
        void setDefaultPatches();

        //- Function to stitch the triangles by removing duplicate points.
        //  Returns true if any points merged
        bool stitchTriangles
        (
            const scalar tol = SMALL,
            const bool verbose = false
        );

        //- Read in OpenFOAM format
        bool readNative(Istream& is);

        //- Write in OpenFOAM format
        void writeNative(Ostream& os) const;

        //- Read in STL format
        bool readSTL(const fileName& filename, bool forceBinary=false);

        //- Generic read routine for given format type.
        //  If the format type is "", uses the file extension.
        bool read
        (
            const fileName& filename,
            const word& fileType,
            const bool check = true
        );

        //- Write STL ASCII format.
        //  Each region becomes a 'solid' 'endsolid' block.
        void writeSTLASCII(const fileName& filename, const bool sort) const;

        //- Write STL BINARY format
        void writeSTLBINARY(const fileName& filename) const;

        //- Write GTS (Gnu Tri Surface library) format.
        void writeGTS(const fileName& filename, const bool sort) const;


    // Static Private Functions

        //- Convert faces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const faceList& faces,
            const label defaultRegion = 0
        );

        //- Convert triFaces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const triFaceList& faces,
            const label defaultRegion = 0
        );

        //- Return a new surface using specified pointMap and faceMap
        //
        //  \param[in] pointMap from subsetMeshMap
        //  \param[in] faceMap from subsetMeshMap
        triSurface subsetMeshImpl
        (
            const labelList& pointMap,
            const labelList& faceMap
        ) const;


protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(MeshReference::points());
        }

        //- Non-const access to the faces
        List<labelledTri>& storedFaces()
        {
            return static_cast<List<labelledTri>&>(*this);
        }


public:

    // Public Typedefs

        //- Placeholder only, but do not remove - it is needed for GeoMesh
        typedef bool BoundaryMesh;

        //- The face type (same as the underlying PrimitivePatch)
        typedef labelledTri face_type;

        //- The point type (same as the underlying PrimitivePatch)
        typedef point point_type;


    //- Runtime type information
    ClassName("triSurface");


    // Static

        //- Name of triSurface directory to use.
        static fileName triSurfInstance(const Time&);

        //- Known readable file-types, including via friends or proxies
        static wordHashSet readTypes();

        //- Known writable file-types, including via friends or proxies
        static wordHashSet writeTypes();

        //- Can we read this file format?
        static bool canReadType(const word& fileType, bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& fileType, bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);


    // IO helpers

        //- Return fileName.
        //  If fileName is relative gets treated as local to IOobject.
        static fileName relativeFilePath
        (
            const IOobject& io,
            const fileName& f,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Return fileName to load IOobject from.
        //  Fatal if the file does not exist
        static fileName checkFile
        (
            const IOobject& io,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Return fileName to load IOobject from.
        //  Supports optional override of fileName with "file" entry
        //  Fatal if the file does not exist
        static fileName checkFile
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Use IOobject information to resolve file to load from,
        //- or empty if the file does not exist.
        static fileName findFile
        (
            const IOobject& io,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Use IOobject information to resolve file to load from,
        //- or empty if the file does not exist.
        //  Supports optional override of fileName with "file" entry
        static fileName findFile
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );


    // Constructors

        //- Default construct
        triSurface();

        //- Copy construct
        triSurface(const triSurface& surf);

        //- Move construct
        triSurface(triSurface&& surf);

        //- Construct from triangles, patches, points.
        triSurface
        (
            const List<labelledTri>& triangles,
            const geometricSurfacePatchList& patches,
            const pointField& pts
        );

        //- Construct from triangles, patches, points. Reuse storage.
        triSurface
        (
            List<labelledTri>& triangles,
            const geometricSurfacePatchList& patches,
            pointField& pts,
            const bool reuse
        );

        //- Construct from triangles, points.
        //- Set patch names to default.
        triSurface
        (
            const List<labelledTri>& triangles,
            const pointField& pts
        );

        //- Construct from triangles, points.
        //- Set region to 0 and default patchName.
        triSurface
        (
            const triFaceList& triangles,
            const pointField& pts
        );

        //- Construct from file name (uses extension to determine type).
        //  Optional (positive, non-zero) point scaling is possible.
        explicit triSurface
        (
            const fileName& name,
            const scalar scaleFactor = -1
        );

        //- Construct from file name with given format type.
        //  If the format type is "", uses the file extension.
        triSurface
        (
            const fileName& name,
            const word& fileType,
            const scalar scaleFactor = -1
        );

        //- Construct from Istream
        explicit triSurface(Istream& is);

        //- Construct from objectRegistry by reading an ".ftr" file
        explicit triSurface(const Time& d);

        //- Read construct using IO to find the file location.
        //  Dictionary may contain the following entries:
        //  - \c file = alternative file name (default is dictionary name)
        //  - \c fileType = file format (default is from file extension)
        //  - \c scale (eg, 0.001: mm to m)
        //  .
        triSurface
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );


    // Selectors

        //- Read construct from filename with given file type
        static autoPtr<triSurface> New
        (
            const fileName& name,
            const word& fileType
        );

        //- Read construct from filename (file type implicit from extension)
        static autoPtr<triSurface> New(const fileName& name);


    //- Destructor
    virtual ~triSurface();


    // Member Functions

        void clearOut();

        void clearTopology();

        void clearPatchMeshAddr();

        void swap(triSurface& surf);


    // Access

        const geometricSurfacePatchList& patches() const noexcept
        {
            return patches_;
        }

        geometricSurfacePatchList& patches() noexcept
        {
            return patches_;
        }

        //- Return const access to the faces
        const List<labelledTri>& surfFaces() const noexcept
        {
            return static_cast<const List<labelledTri>&>(*this);
        }

        //- Return edge-face addressing sorted (for edges with more than
        //  2 faces) according to the angle around the edge.
        //  Orientation is anticlockwise looking from
        //  edge.vec(localPoints())
        const labelListList& sortedEdgeFaces() const;

        //- If 2 face neighbours: label of face where ordering of edge
        //  is consistent with righthand walk.
        //  If 1 neighbour: label of only face.
        //  If >2 neighbours: undetermined.
        const labelList& edgeOwner() const;


        //- Face area vectors (normals)
        const vectorField& Sf() const
        {
            return MeshReference::faceAreas();
        }

        //- Face area magnitudes
        const scalarField& magSf() const
        {
            return MeshReference::magFaceAreas();
        }

        //- Face centres
        const vectorField& Cf() const
        {
            return MeshReference::faceCentres();
        }


    // Interoperability with other surface mesh classes

        //- Sort faces according to zoneIds
        //  Returns a surfZoneList and sets faceMap to index within faces()
        //  (i.e. map from original,unsorted to sorted)
        List<surfZone> sortedZones(labelList& faceMap) const;

        //- Create a list of faces from the triFaces
        void triFaceFaces(List<face>& plainFaceList) const;


    // Edit

        //- Move points
        virtual void movePoints(const pointField& pts);

        //- Swap points. Similar to movePoints, but returns the old points
        virtual void swapPoints(pointField& pts);

        //- Scale points. A non-positive factor is ignored.
        virtual void scalePoints(const scalar scaleFactor);

        //- Check/remove duplicate/degenerate triangles
        void checkTriangles(const bool verbose);

        //- Check triply (or more) connected edges.
        void checkEdges(const bool verbose);

        //- Remove non-valid triangles
        void cleanup(const bool verbose);

        //- Remove unused points and renumber faces in local visit order
        //
        //  \param[out] pointMap from new to old points (optional)
        void compactPoints
        (
            labelList& pointMap = const_cast<labelList&>(labelList::null())
        );

        //- Fill faceZone with currentZone for every face reachable
        //  from facei without crossing edge marked in borderEdge.
        //  Note: faceZone has to be sized nFaces before calling this fun.
        void markZone
        (
            const boolList& borderEdge,
            const label facei,
            const label currentZone,
            labelList& faceZone
        ) const;

        //- (size and) fills faceZone with zone of face. Zone is area
        //  reachable by edge crossing without crossing borderEdge
        //  (bool for every edge in surface). Returns number of zones.
        label markZones
        (
            const boolList& borderEdge,
            labelList& faceZone
        ) const;

        //- Create mappings for a sub-surface
        //
        //  \param[in] include the faces to select
        //  \param[out] pointMap from new to old localPoints
        //  \param[out] faceMap from new to old localFaces
        template<class BoolListType>
        void subsetMeshMap
        (
            const BoolListType& include,
            labelList& pointMap,
            labelList& faceMap
        ) const
        {
            PatchTools::subsetMap(*this, include, pointMap, faceMap);
        }

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        //  \param[out] pointMap from new to old localPoints
        //  \param[out] faceMap from new to old localFaces
        triSurface subsetMesh
        (
            const UList<bool>& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        //  \param[out] pointMap from subsetMeshMap
        //  \param[out] faceMap from subsetMeshMap
        triSurface subsetMesh
        (
            const bitSet& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        triSurface subsetMesh(const UList<bool>& include) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        triSurface subsetMesh(const bitSet& include) const;

        //- Return a new surface subsetted on the selected patch names
        //
        //  \param[in] includeNames surface patch names to include
        //  \param[in] excludeNames surface patch names to exclude
        //
        //  \see Foam::stringListOps::findMatching for details about matching
        triSurface subsetMesh
        (
            const wordRes& includeNames,
            const wordRes& excludeNames = wordRes()
        ) const;

        //- Swap the list of faces being addressed
        void swapFaces(List<labelledTri>& faceLst);

        //- Alter contents by transferring (triangles, points) components.
        //  Patch information is small and therefore just copied.
        void transfer(triSurface& surf);

        //- Alter contents by transferring (triangles, points) components.
        //  Patch information is small and therefore just copied.
        void transfer(MeshedSurface<labelledTri>& surf);


    // Write

        //- Write to Ostream in simple OpenFOAM format
        void write(Ostream& os) const;

        //- Generic write routine (uses extension to determine type).
        //  The sort option may not have an effect.
        void write(const fileName&, const bool sortByRegion = false) const;

        //- Generic write routine for given format type.
        //  If the format type is "", uses the file extension.
        //  The sort option may not have an effect.
        void write
        (
            const fileName& filename,
            const word& fileType,
            const bool sortByRegion = false
        ) const;

        //- Write to database
        void write(const Time& d) const;

        //- Write some statistics
        void writeStats(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const triSurface& surf);

        //- Move assignment
        void operator=(triSurface&& surf);

        //- Move assignment
        void operator=(MeshedSurface<labelledTri>&& surf);


    // IOstream Operators

        friend Istream& operator>>(Istream& is, triSurface& s);
        friend Ostream& operator<<(Ostream& os, const triSurface& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
