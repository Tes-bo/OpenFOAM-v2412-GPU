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
    Foam::surfMesh

Description
    A surface mesh consisting of general polygon faces that has
    IO capabilities and a registry for storing fields.

Note
    In many places the Foam::polySurface class may be more flexible.

SourceFiles
    surfMesh.C
    surfMeshClear.C
    surfMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfMesh_H
#define Foam_surfMesh_H

#include "surfaceRegistry.H"
#include "MeshedSurfaceIOAllocator.H"
#include "primitivePatch.H"
#include "surfZoneIOList.H"
#include "surfFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dimensionSet;
class meshedSurf;
template<class Face> class MeshedSurface;

/*---------------------------------------------------------------------------*\
                          Class surfMesh Declaration
\*---------------------------------------------------------------------------*/

class surfMesh
:
    public  surfaceRegistry,
    private Detail::MeshedSurfaceIOAllocator,
    public  PrimitivePatch<::Foam::UList<face>, const pointField&>
{
public:

    // Public Data Types

        //- Enumeration defining the state of the mesh after a read update.
        //  Used for post-processing applications, where the mesh
        //  needs to update based on the files written in time directories
        enum readUpdateState
        {
            UNCHANGED,
            POINTS_MOVED,
            TOPO_CHANGE,
            TOPO_PATCH_CHANGE
        };


private:

    // Private Typedefs

        //- Memory and IO management
        typedef
            Detail::MeshedSurfaceIOAllocator
            Allocator;

        //- Internal mesh storage type
        typedef
            PrimitivePatch<::Foam::UList<face>, const pointField&>
            MeshReference;


    // Private Data

        //- The surface zones
        surfZoneIOList surfZones_;


protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return Allocator::storedIOPoints();
        }

        //- Non-const access to the faces
        faceList& storedFaces()
        {
            return Allocator::storedIOFaces();
        }

        //- Non-const access to the surface zones
        surfZoneList& storedZones()
        {
            return surfZones_;
        }

        //- Update point/face references
        void updateRefs();


public:

    // Public Typedefs

    //- Typedef required for GeoMesh
    typedef surfMesh Mesh;

    //- Placeholder only, but do not remove - it is needed for GeoMesh
    typedef bool BoundaryMesh;


    //- Declare type-name, virtual type (with debug switch)
    TypeName("surfMesh");

    //- Return the mesh sub-directory name (normally "surfMesh")
    static word meshSubDir;


    // Generated Methods

        //- No copy construct
        surfMesh(const surfMesh&) = delete;

        //- No copy assignment
        void operator=(const surfMesh&) = delete;


    // Constructors

        //- Read construct from IOobject.
        //  Writing = NO_WRITE
        explicit surfMesh(const IOobject& io);

        //- Read construct from IOobject, with alternative surface name
        //  Writing = NO_WRITE
        surfMesh(const IOobject& io, const word& surfName);

        //- Construct empty with specified name on the given registry.
        surfMesh(const word& surfName, const objectRegistry& obr);

        //- Copy construct from MeshedSurface<face>
        surfMesh
        (
            const IOobject& io,
            const MeshedSurface<face>& surf,
            const word& surfName = word::null
        );

        //- Move construct from MeshedSurface<face>
        surfMesh
        (
            const IOobject& io,
            MeshedSurface<face>&& surf,
            const word& surfName = word::null
        );


    //- Destructor
    virtual ~surfMesh();


    // Member Functions

    // Database

        //- Return the local mesh directory (dbDir()/meshSubDir)
        fileName meshDir() const;

        //- Return the current instance directory for points
        //  Used when constructing geometric mesh data dependent on points
        const fileName& pointsInstance() const;

        //- Return the current instance directory for faces
        const fileName& facesInstance() const;

        //- Set the instance for mesh files
        void setInstance
        (
            const fileName& inst,
            IOobjectOption::writeOption wOpt = IOobjectOption::AUTO_WRITE
        );

        //- Adjust the write option for all components
        void setWriteOption(IOobjectOption::writeOption wOpt);


    // Access

        //- Return the number of raw points
        virtual label nPoints() const;

        //- Return the number of raw faces
        virtual label nFaces() const;

        //- Return number of faces
        virtual label size() const
        {
            return nFaces();
        }


        //- Return points
        virtual const pointField& points() const;

        //- Return faces
        virtual const faceList& faces() const;

        //- Return surface zones
        virtual const surfZoneList& surfZones() const
        {
            return surfZones_;
        }

        //- Return face area vectors (normals)
        const vectorField& Sf() const
        {
            return MeshReference::faceAreas();
        }

        //- Return face area magnitudes
        const scalarField& magSf() const
        {
            return MeshReference::magFaceAreas();
        }

        //- Face centres
        const vectorField& Cf() const
        {
            return MeshReference::faceCentres();
        }


    // Zones

        //- Add surface zones, optionally validating the zone coverage
        void addZones
        (
            const surfZoneList& zones,
            bool validate = true
        );

        //- Remove surface zones
        void removeZones();

        //- Check the surface zone definitions
        void checkZones(const bool verbose=true);


    // Modification

        //- Update with new contents
        void copySurface
        (
            const pointField& points,
            const faceList& faces,
            bool validate=false
        );

        //- Update with new contents
        void copySurface
        (
            const meshedSurf& surf,
            bool validate=false
        );

        //- Update with new contents
        void copySurface
        (
            const MeshedSurface<face>& surf,
            bool validate=false
        );

        //- Transfer the contents of the argument and annul the argument
        //  Optionally validate the zone coverage.
        void transfer
        (
            MeshedSurface<face>& surf,
            bool validate=false
        );

        //- Update mesh based on the files saved in time directories
        virtual readUpdateState readUpdate();


    // Fields

        //- Create/store named zero field as face or point data
        //- (template parameter).
        //
        //  - Default is face-data (surfGeoMesh as template).
        //  - For point-data use surfPointGeoMesh as template.
        //
        //  \return reference to the field
        template<class Type, class GeoMeshType = surfGeoMesh>
        DimensionedField<Type, GeoMeshType>&
        newField
        (
            const word& fieldName,
            const dimensionSet& dims
        );

        //- Copy/store named field as face or point data (template parameter).
        //
        //  - Default is face-data (surfGeoMesh as template).
        //  - For point-data use surfPointGeoMesh as template.
        //
        //  \return reference to the field
        template<class Type, class GeoMeshType = surfGeoMesh>
        DimensionedField<Type, GeoMeshType>&
        storeField
        (
            const word& fieldName,
            const dimensionSet& dims,
            const Field<Type>& values
        );

        //- Move/store named field as face or point data (template parameter).
        //
        //  - Default is face-data (surfGeoMesh as template).
        //  - For point-data use surfPointGeoMesh as template.
        //
        //  \return reference to the field
        template<class Type, class GeoMeshType = surfGeoMesh>
        DimensionedField<Type, GeoMeshType>&
        storeField
        (
            const word& fieldName,
            const dimensionSet& dims,
            Field<Type>&& values
        );


    // Writing

        //- Avoid masking the normal objectRegistry write
        using surfaceRegistry::write;

        //- Write all components using given format, version and compression
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;


        //- Write to file, choosing writer based on its extension.
        //  Uses MeshedSurfaceProxy for writing.
        void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const;

        //- Write to file, choosing writer for given fileType.
        //  Uses MeshedSurfaceProxy for writing.
        void write
        (
            const fileName& name,
            const word& fileType,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const;


    // Storage Management

        //- Release the geometry and return as a MeshedSurface<face>.
        autoPtr<MeshedSurface<face>> releaseGeom();

        //- Clear geometry
        void clearGeom();

        //- Clear addressing
        void clearAddressing();

        //- Clear all geometry and addressing unnecessary for CFD
        void clearOut();

        //- Clear primitive data (points, faces and cells)
        void clearPrimitives();

        //- Clear stored fields
        void clearFields();

        //- Remove all files from mesh instance
        void removeFiles(const fileName& instanceDir) const;

        //- Remove all files from mesh instance()
        void removeFiles() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
