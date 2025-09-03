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
    Foam::MeshedSurface

Description
    A surface geometry mesh with zone information, not to be confused with
    the similarly named surfaceMesh, which actually refers to the cell faces
    of a volume mesh.

    A MeshedSurface can have zero or more surface zones (roughly equivalent
    to faceZones for a polyMesh). If surface zones are defined, they must
    be contiguous and cover all of the faces.

    The MeshedSurface is intended for surfaces from a variety of sources.
    - A set of points and faces without any surface zone information.
    - A set of points and faces with randomly ordered zone information.
      This could arise, for example, from reading external file formats
      such as STL, etc.

SourceFiles
    MeshedSurface.C
    MeshedSurfaceCore.C
    MeshedSurfaceIO.C
    MeshedSurfaceNew.C
    MeshedSurfaceZones.C
    MeshedSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_MeshedSurface_H
#define Foam_MeshedSurface_H

#include "primitivePatch.H"
#include "PatchTools.H"
#include "face.H"
#include "labelledTri.H"
#include "bitSet.H"
#include "HashSet.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class surfMesh;
class polyBoundaryMesh;

template<class Face> class MeshedSurface;
template<class Face> class MeshedSurfaceProxy;
template<class Face> class UnsortedMeshedSurface;

template<class Face>
Istream& operator>>(Istream&, MeshedSurface<Face>&);
template<class Face>
Ostream& operator<<(Ostream&, const MeshedSurface<Face>&);

/*---------------------------------------------------------------------------*\
                      Class MeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurface
:
    public PrimitivePatch<::Foam::List<Face>, pointField>,
    public fileFormats::surfaceFormatsCore
{
    // Friends, regardless of face representations
    template<class Face2> friend class MeshedSurface;
    template<class Face2> friend class UnsortedMeshedSurface;

    // Friendship with surfMesh is needed for transferring
    friend class surfMesh;


private:

    // Private Typedefs (convenience)

        //- Internal mesh storage type
        typedef PrimitivePatch<::Foam::List<Face>, pointField>
            MeshReference;

        typedef UnsortedMeshedSurface<Face>  FriendType;
        typedef MeshedSurfaceProxy<Face>     ProxyType;


    // Private Data

        //- Face ids.
        //  If these exist, they are typically arise from reading a mesh
        //  format from another CAE software (eg, NASTRAN, STARCD, ...)
        labelList faceIds_;

        //- Zone information
        // (face ordering nFaces/startFace only used during reading/writing)
        surfZoneList zones_;


    // Private Member Functions

        //- Read/construct from Istream
        Istream& read(Istream& is);

        //- Write to Ostream
        Ostream& write(Ostream& os) const;

        //- Return a new surface using specified pointMap and faceMap
        //
        //  \param[in] pointMap from subsetMeshMap
        //  \param[in] faceMap from subsetMeshMap
        MeshedSurface subsetMeshImpl
        (
            const labelList& pointMap,
            const labelList& faceMap
        ) const;


protected:

    // Protected Member Functions

        //- Transfer points/zones from 'face' to other other shapes.
        //  Eg, transcribe face to triFace, or face -> labelledTri, including
        //  any addZonesToFaces adjustment.
        //  No general form, only specializations.
        void transcribe(MeshedSurface<face>& surf);

        //- Sanity check/resizing on zones.
        //  Adjust zones so that they cover the number of faces
        //  The last zone will be extended as needed
        void checkZones(const bool verbose = true);

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(MeshReference::points());
        }

        //- Non-const access to the faces
        List<Face>& storedFaces()
        {
            return static_cast<List<Face>&>(*this);
        }

        //- Non-const access to face ids
        labelList& storedFaceIds()
        {
            return faceIds_;
        }

        //- Non-const access to the zones
        surfZoneList& storedZones()
        {
            return zones_;
        }

        //- Sort faces by zones and store sorted faces
        void sortFacesAndStore
        (
            DynamicList<Face>& unsortedFaces,
            DynamicList<label>& zoneIds,
            DynamicList<label>& elemIds,
            bool sorted
        );

        //- Set new zones from faceMap
        virtual void remapFaces(const labelUList& faceMapNewToOld);


public:

    // Public Typedefs

        //- The face type (same as the underlying PrimitivePatch)
        typedef Face face_type;

        //- The point type (same as the underlying PrimitivePatch)
        typedef point point_type;


    //- Declare type-name (with debug switch)
    ClassName("MeshedSurface");


    // Static Functions

        //- Known readable file-types, without friends or proxies
        static wordHashSet readTypes();

        //- Known writable file-types, without friends or proxies
        static wordHashSet writeTypes();

        //- Can we read this file format? Also checks friend types.
        static bool canReadType(const word& fileType, bool verbose=false);

        //- Can we write this file format? Also checks proxy types.
        static bool canWriteType(const word& fileType, bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);


    // Constructors

        //- Default construct, an empty surface
        MeshedSurface();

        //- Copy construct
        MeshedSurface(const MeshedSurface& surf);

        //- Copy construct from an UnsortedMeshedSurface
        MeshedSurface(const UnsortedMeshedSurface<Face>& surf);

        //- Move construct
        MeshedSurface(MeshedSurface&& surf);

        //- Move construct from an UnsortedMeshedSurface
        MeshedSurface(UnsortedMeshedSurface<Face>&& surf);

        //- Copy construct from components (points, faces, zones).
        MeshedSurface
        (
            const pointField& pointLst,
            const UList<Face>& faceLst,
            const UList<surfZone>& zoneLst
        );

        //- Move construct from components (points, faces).
        //  Zone information is fairly lightweight and is copied.
        MeshedSurface
        (
            pointField&& pointLst,
            List<Face>&& faceLst,
            const UList<surfZone>& zoneLst
        );

        //- Copy construct from components (points, faces).
        //  Use zone information if available
        MeshedSurface
        (
            const pointField& pointLst,
            const UList<Face>& faceLst,
            const labelUList& zoneSizes = labelUList(),
            const UList<word>& zoneNames = UList<word>()
        );

        //- Move construct from components (points, faces).
        //  Use zone information if available
        MeshedSurface
        (
            pointField&& pointLst,
            List<Face>&& faceLst,
            const labelUList& zoneSizes = labelUList(),
            const UList<word>& zoneNames = UList<word>()
        );

        //- Construct from a boundary mesh with local points/faces
        MeshedSurface
        (
            const polyBoundaryMesh& bMesh,
            const bool globalPoints = false
        );

        //- Construct from a surfMesh
        explicit MeshedSurface(const surfMesh& mesh);

        //- Construct from file name (uses extension to determine type)
        explicit MeshedSurface(const fileName& name);

        //- Construct from file name and given file type
        //  If the format type is "", uses the file extension.
        explicit MeshedSurface(const fileName& name, const word& fileType);

        //- Construct from Istream
        explicit MeshedSurface(Istream& is);

        //- Construct from database (as surfMesh) with default name
        explicit MeshedSurface(const Time& runTime);

        //- Construct from database (as surfMesh) with given surface name
        MeshedSurface(const Time& runTime, const word& surfName);

        //- Read construct using IO to find the file location.
        //  Dictionary may contain the following entries:
        //  - \c file = alternative file name (default is dictionary name)
        //  - \c fileType = file format (default is from file extension)
        //  - \c scale (eg, 0.001: mm to m)
        //  .
        MeshedSurface
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            MeshedSurface,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Read construct from filename with given file type
        //
        //  \note Use mandatory=false if support for the file type
        //  is optional (the file still needs to exist!).
        static autoPtr<MeshedSurface> New
        (
            const fileName& name,
            const word& fileType,
            bool mandatory = true
        );

        //- Read construct from filename (file type implicit from extension)
        static autoPtr<MeshedSurface> New(const fileName& name);


    //- Destructor
    virtual ~MeshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& name,
                const MeshedSurface<Face>& surf,
                IOstreamOption streamOpt,
                const dictionary& options
            ),
            (name, surf, streamOpt, options)
        );

        //- Write to file, selecting writer based on its extension
        static void write
        (
            const fileName& name,
            const MeshedSurface<Face>& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );

        //- Write to file, selecting writer based on the given extension
        static void write
        (
            const fileName& name,
            const word& fileType,
            const MeshedSurface<Face>& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

    // Access

        //- The surface size is the number of faces
        label size() const
        {
            return MeshReference::size();
        }

        //- Return const access to the faces
        const List<Face>& surfFaces() const
        {
            return static_cast<const List<Face>&>(*this);
        }

        //- Return const access to faces ids
        //  If these exist, they are typically arise from reading a mesh
        //  format from another CAE software (eg, NASTRAN, STARCD, ...)
        const labelList& faceIds() const
        {
            return faceIds_;
        }

        //- Const access to the surface zones.
        //  If zones are defined, they must be contiguous and cover the
        //  entire surface
        const surfZoneList& surfZones() const
        {
            return zones_;
        }

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


    // Edit

        //- Clear all storage
        virtual void clear();


        //- Add surface zones
        virtual void addZones
        (
            const UList<surfZone>&,
            const bool cullEmpty=false
        );

        //- Add surface zones
        virtual void addZones
        (
            const labelUList& sizes,
            const UList<word>& names,
            const bool cullEmpty=false
        );

        //- Add surface zones
        virtual void addZones
        (
            const labelUList& sizes,
            const bool cullEmpty=false
        );

        //- Propagate zone information on face regions.
        //  Normally a no-op, only used by the labelledTri specialization.
        //  Specializations return true, others return false.
        bool addZonesToFaces();


        //- Remove surface zones
        virtual void removeZones();


        //- Move points
        virtual void movePoints(const pointField& newPoints);

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar scaleFactor);

        //- Remove invalid faces
        virtual void cleanup(const bool verbose);

        //- Remove unused points and renumber faces in local visit order
        //
        //  \param[out] pointMap from new to old points (optional)
        virtual void compactPoints
        (
            labelList& pointMap = const_cast<labelList&>(labelList::null())
        );

        virtual bool stitchFaces
        (
            const scalar tol=SMALL,
            const bool verbose=false
        );

        virtual bool checkFaces
        (
            const bool verbose=false
        );

        //- Count number of triangles.
        virtual label nTriangles() const;

        //- Count number of triangles, returning a face map of original ids.
        //  The faceMap is zero-sized when no triangulation would be needed.
        virtual label nTriangles(labelList& faceMap) const;

        //- Triangulate in-place, returning the number of triangles added.
        virtual label triangulate();

        //- Triangulate in-place, returning the number of triangles added
        //  and setting a map of original face Ids.
        //  The faceMap is zero-sized when no triangulation was done.
        virtual label triangulate(labelList& faceMap);

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
        MeshedSurface subsetMesh
        (
            const UList<bool>& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        //  \param[out] pointMap from new to old localPoints
        //  \param[out] faceMap from new to old localFaces
        MeshedSurface subsetMesh
        (
            const bitSet& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        MeshedSurface subsetMesh(const UList<bool>& include) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        MeshedSurface subsetMesh(const bitSet& include) const;

        //- Return a new surface subsetted on the selected zone names
        //
        //  \param[in] includeNames surface zone names to include
        //  \param[in] excludeNames surface zone names to exclude
        //
        //  \see Foam::stringListOps::findMatching for details about matching
        MeshedSurface subsetMesh
        (
            const wordRes& includeNames,
            const wordRes& excludeNames = wordRes()
        ) const;

        //- Swap contents
        void swap(MeshedSurface<Face>& surf);

        //- Transfer the components
        void transfer(pointField& pointLst, List<Face>& faceLst);

        //- Transfer the contents of the argument and annul the argument
        void transfer(MeshedSurface<Face>& surf);

        //- Transfer the contents of the argument and annul the argument
        void transfer(UnsortedMeshedSurface<Face>& surf);

        //- Release (clear) geometry and return for reuse
        autoPtr<MeshedSurface<Face>> releaseGeom();

        //- Swap the stored faces. Use with caution
        void swapFaces(List<Face>& faces);

        //- Swap the stored points
        void swapPoints(pointField& points);


    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName& name, const word& fileType);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);


    // Write

        void writeStats(Ostream& os) const;

        //- Generic write routine. Chooses writer based on extension.
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, *this, streamOpt, options);
        }

        //- Generic write routine for given format type.
        //  If the format type is "", uses the file extension.
        virtual void write
        (
            const fileName& name,
            const word& fileType,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, fileType, *this, streamOpt, options);
        }


        //- Write to database
        void write
        (
            const Time& runTime,
            const word& surfName = word::null
        ) const;


    // Member Operators

        //- Copy assignment
        void operator=(const MeshedSurface<Face>& surf);

        //- Move assignment
        void operator=(MeshedSurface<Face>&& surf);

        //- Conversion operator to MeshedSurfaceProxy
        operator MeshedSurfaceProxy<Face>() const;


    // IOstream Operators

        //- Read MeshedSurface from Istream.
        //  Avoid using to read/write file content (fragile).
        friend Istream& operator>> <Face>
        (
            Istream& is,
            MeshedSurface<Face>& surf
        );


        //- Write MeshedSurface to Ostream.
        //  Avoid using to read/write file content (fragile).
        friend Ostream& operator<< <Face>
        (
            Ostream& os,
            const MeshedSurface<Face>& surf
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for labelledTri.
template<>
bool MeshedSurface<labelledTri>::addZonesToFaces();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
