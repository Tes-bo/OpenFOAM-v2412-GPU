/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::UnsortedMeshedSurface

Description
    A surface geometry mesh, in which the surface zone information is
    conveyed by the 'zoneId' associated with each face.

    This form of surface description is particularly useful for reading in
    surface meshes from third-party formats (eg, obj, stl, gts, etc.). It
    can also be particularly useful for situations in which the surface
    many be adjusted in an arbitrary manner without worrying about needed
    to adjust the zone information (eg, surface refinement).

See also
    The Foam::MeshedSurface - which is organized as a surface mesh, but
    with independent zone information.

SourceFiles
    UnsortedMeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UnsortedMeshedSurface_H
#define Foam_UnsortedMeshedSurface_H

#include "MeshedSurface.H"
#include "surfZoneIdentifierList.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "HashSet.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

template<class Face> class MeshedSurface;
template<class Face> class MeshedSurfaceProxy;
template<class Face> class UnsortedMeshedSurface;

template<class Face>
Istream& operator>>(Istream&, UnsortedMeshedSurface<Face>&);
template<class Face>
Ostream& operator<<(Ostream&, const UnsortedMeshedSurface<Face>&);

/*---------------------------------------------------------------------------*\
                   Class UnsortedMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class UnsortedMeshedSurface
:
    public MeshedSurface<Face>
{
    // Friends, regardless of face representations
    template<class Face2> friend class MeshedSurface;
    template<class Face2> friend class UnsortedMeshedSurface;

    // Friendship with surfMesh is needed for transferring
    friend class surfMesh;

private:

    // Private Typedefs (convenience)

        typedef MeshedSurface<Face>       MeshReference;
        typedef MeshedSurface<Face>       FriendType;
        typedef MeshedSurfaceProxy<Face>  ProxyType;


    // Private Data

        //- The zone Id associated with each face
        labelList zoneIds_;

        //- Zone information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfZoneIdentifier> zoneToc_;


    // Private Member Functions

        //- Disable resize with value
        void resize(const label, const Face&) = delete;

        //- Disable setSize with value
        void setSize(const label, const Face&) = delete;

        //- Read/construct from Istream
        bool readIstream(Istream& is);

        //- Write to Ostream
        void writeOstream(Ostream& os) const;


        //- Return a new surface using specified pointMap and faceMap
        //
        //  \param[in] pointMap from subsetMeshMap
        //  \param[in] faceMap from subsetMeshMap
        UnsortedMeshedSurface subsetMeshImpl
        (
            const labelList& pointMap,
            const labelList& faceMap
        ) const;


protected:

    // Protected Member Functions

        //- Return non-const access to the zone Ids
        labelList& storedZoneIds()
        {
            return zoneIds_;
        }

        //- Return non-const access to the zone table-of-contents
        List<surfZoneIdentifier>& storedZoneToc()
        {
            return zoneToc_;
        }

        //- Set new zones from faceMap
        virtual void remapFaces(const labelUList& faceMapNewToOld);


public:

    // Public Typedefs

        //- The face type (same as the underlying PrimitivePatch)
        typedef Face face_type;

        //- The point type (same as the underlying PrimitivePatch)
        typedef point point_type;


    //- Declare type-name, virtual type (with debug switch)
    TypeName("UnsortedMeshedSurface");


    // Static Functions

        //- Known readable file-types, without friends or proxies
        static wordHashSet readTypes();

        //- Known writable file-types, without friends or proxies
        static wordHashSet writeTypes();

        //- Can we read this file format? Also checks friend types.
        static bool canReadType(const word& fileType, bool verbose=false);

        //- Can we write this file format? Also checks friend types.
        static bool canWriteType(const word& fileType, bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);


    // Constructors

        //- Default construct
        UnsortedMeshedSurface();

        //- Copy construct
        UnsortedMeshedSurface(const UnsortedMeshedSurface<Face>& surf);

        //- Construct from a MeshedSurface
        UnsortedMeshedSurface(const MeshedSurface<Face>& surf);

        //- Move construct from a UnsortedMeshedSurface
        UnsortedMeshedSurface(UnsortedMeshedSurface<Face>&& surf);

        //- Move construct from MeshedSurface
        UnsortedMeshedSurface(MeshedSurface<Face>&& surf);

        //- Move construct from components (points, faces, zone ids).
        //  Zone information is fairly lightweight and is copied.
        UnsortedMeshedSurface
        (
            pointField&& pointLst,
            List<Face>&& faceLst,
            List<label>&& zoneIds,
            UList<surfZoneIdentifier>& tocInfo
        );

        //- Construct from file name (uses extension to determine type)
        explicit UnsortedMeshedSurface(const fileName& name);

        //- Construct from file name with given format type.
        //  If the format type is "", uses the file extension.
        UnsortedMeshedSurface(const fileName& name, const word& fileType);

        //- Construct from Istream
        explicit UnsortedMeshedSurface(Istream& is);

        //- Construct from database (as surfMesh) with default name
        explicit UnsortedMeshedSurface(const Time& runTime);

        //- Construct from database (as surfMesh) with given surface name
        UnsortedMeshedSurface(const Time& runTime, const word& surfName);

        //- Read construct using IO to find the file location.
        //  Dictionary may contain the following entries:
        //  - \c file = alternative file name (default is dictionary name)
        //  - \c fileType = file format (default is from file extension)
        //  - \c scale (eg, 0.001: mm to m)
        //  .
        UnsortedMeshedSurface
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            UnsortedMeshedSurface,
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
        static autoPtr<UnsortedMeshedSurface> New
        (
            const fileName& name,
            const word& fileType,
            bool mandatory = true
        );

        //- Read construct from filename (implicit extension)
        static autoPtr<UnsortedMeshedSurface> New(const fileName& name);


    //- Destructor
    virtual ~UnsortedMeshedSurface() = default;


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& name,
                const UnsortedMeshedSurface<Face>& surf,
                IOstreamOption streamOpt,
                const dictionary& options
            ),
            (name, surf, streamOpt, options)
        );

        //- Write to file, select based on its extension
        static void write
        (
            const fileName& name,
            const UnsortedMeshedSurface<Face>& surf,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );

        //- Write to file with given format type.
        //  If the format type is "", uses the file extension.
        static void write
        (
            const fileName& name,
            const word& fileType,
            const UnsortedMeshedSurface<Face>& surf,
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

        //- Reset size of face and zone list
        void setSize(const label);

        //- Return const access to the zone ids
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }

        //- Return const access to the zone table-of-contents
        const List<surfZoneIdentifier>& zoneToc() const
        {
            return zoneToc_;
        }

        //- Sort faces according to zoneIds
        //  Returns a surfZoneList and sets faceMap to index within faces()
        //  (i.e. map from original,unsorted to sorted)
        surfZoneList sortedZones(labelList& faceMap) const;

        //- Set zones to 0 and set a single zone
        void setOneZone();

        //- Set zone ids and zones
        void setZones(const surfZoneList& zoneLst);

        //- Set zone ids and zones
        void setZones(const labelUList& sizes, const UList<word>& names);

        //- Set zone ids and zones with default names
        void setZones(const labelUList& sizes);


    // Edit

        //- Clear all storage
        virtual void clear();

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
        UnsortedMeshedSurface subsetMesh
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
        UnsortedMeshedSurface subsetMesh
        (
            const bitSet& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        UnsortedMeshedSurface subsetMesh(const UList<bool>& include) const;

        //- Return a new surface subsetted on the selected faces.
        //
        //  \param[in] include the faces to select
        UnsortedMeshedSurface subsetMesh(const bitSet& include) const;


        //- Swap contents - disabled
        void swap(MeshedSurface<Face>& surf) = delete;

        //- Swap contents
        void swap(UnsortedMeshedSurface<Face>& surf);

        //- Transfer the contents of the argument and annul the argument
        void transfer(UnsortedMeshedSurface<Face>& surf);

        //- Transfer the contents of the argument and annul the argument
        void transfer(MeshedSurface<Face>& surf);

        //- Release (clear) stored zoneIds and return for reuse
        autoPtr<labelList> releaseZoneIds();


    // Read

        //- Read from file with given format type.
        //  If the format type is "", uses the file extension.
        bool read(const fileName& name, const word& fileType);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);


    // Write

        //- Write to file, choosing writer based on the file extension.
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, *this, streamOpt, options);
        }

        //- Write to file with given format type.
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
            const Time& t,
            const word& surfName = word::null
        ) const;


    // Member Operators

        //- Copy assignment
        void operator=(const UnsortedMeshedSurface<Face>& surf);

        //- Move assignment
        void operator=(UnsortedMeshedSurface<Face>&& surf);

        //- Conversion operator to MeshedSurfaceProxy
        operator MeshedSurfaceProxy<Face>() const;


    // IOstream Operators

        //- Read UnsortedMeshedSurface from Istream.
        //  Avoid using to read/write file content (fragile).
        friend Istream& operator>> <Face>
        (
            Istream& is,
            UnsortedMeshedSurface<Face>& surf
        );

        //- Write UnsortedMeshedSurface to Ostream.
        //  Avoid using to read/write file content (fragile).
        friend Ostream& operator<< <Face>
        (
            Ostream& os,
            const UnsortedMeshedSurface<Face>& surf
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UnsortedMeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
