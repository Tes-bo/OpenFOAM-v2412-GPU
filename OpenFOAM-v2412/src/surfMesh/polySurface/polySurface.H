/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::polySurface

Description
    A surface mesh consisting of general polygon faces and capable of
    holding fields.

SourceFiles
    polySurface.C
    polySurfaceClear.C
    polySurfaceIO.C
    polySurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_polySurface_H
#define Foam_polySurface_H

#include "objectRegistry.H"
#include "primitivePatch.H"
#include "meshedSurf.H"
#include "polySurfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dimensionSet;
class surfZone;
class polySurfaceGeoMesh;
class polySurfacePointGeoMesh;

template<class Face> class MeshedSurface;


/*---------------------------------------------------------------------------*\
                         Class polySurface Declaration
\*---------------------------------------------------------------------------*/

class polySurface
:
    public objectRegistry,
    public PrimitivePatch<::Foam::List<face>, pointField>,
    public meshedSurf
{
public:

    //- Enumeration for the field association
    enum FieldAssociation
    {
        NO_DATA = 0,            //!< No associated data
        FACE_DATA = 1,          //!< Data associated with faces
        POINT_DATA = 2,         //!< Data associated with points
        FACE_OR_POINT_DATA = 3  //!< Data associated with faces or points
    };


private:

    // Private Typedefs

        //- Internal mesh storage type
        typedef
            PrimitivePatch<::Foam::List<face>, pointField>
            MeshReference;


    // Private Data

        //- Per-face zone/region information
        labelList zoneIds_;


    // Private Member Functions

        //- Calculate per-face zone/region information
        void calculateZoneIds(const UList<surfZone>& zones);


protected:

    // Protected Member Functions

        //- Non-const access to points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(MeshReference::points());
        }

        //- Non-const access to the faces
        faceList& storedFaces()
        {
            return
                static_cast<faceList&>
                (
                    static_cast<MeshReference&>(*this)
                );
        }

        //- Const access to the faces
        const faceList& storedFaces() const
        {
            return
                static_cast<const faceList&>
                (
                    static_cast<const MeshReference&>(*this)
                );
        }


public:

    // Public Typedefs

    //- Typedef required for GeoMesh
    typedef polySurface Mesh;

    //- Placeholder only, for GeoMesh
    typedef bool BoundaryMesh;

    //- Name for point fields sub-registry
    static const word pointDataName;


    //- Runtime type information
    TypeName("polySurface");


    // Generated Methods

        //- No copy construct
        polySurface(const polySurface&) = delete;

        //- No copy assignment
        void operator=(const polySurface&) = delete;


    // Constructors

        //- Default construct with NO_READ, NO_WRITE
        //- optionally with a checkIn on the parent registry.
        //  Created without a PointData sub-registry
        explicit polySurface(const IOobject& io, bool doCheckIn = false);

        //- Default construct with specified name on the given registry,
        //- optionally with a checkIn on the parent registry.
        //  Created without a PointData sub-registry
        polySurface
        (
            const word& surfName,
            const objectRegistry& obr,
            bool doCheckIn = false
        );

        //- Copy construct from MeshedSurface<face> contents
        //- with NO_READ, NO_WRITE
        //- optionally with a checkIn on the parent registry.
        //  Created without a PointData sub-registry
        polySurface
        (
            const IOobject& io,
            const MeshedSurface<face>& surf,
            bool doCheckIn = false
        );

        //- Move construct from MeshedSurface<face> contents
        //- with NO_READ, NO_WRITE
        //- optionally with a checkIn on the parent registry.
        //  Created without a PointData sub-registry
        polySurface
        (
            const IOobject& io,
            MeshedSurface<face>&& surf,
            bool doCheckIn = false
        );


    //- Destructor
    virtual ~polySurface();


    // Factory Methods

        //- Get or create (NO_READ, NO_WRITE) named surface on registry
        static polySurface& New(const word& surfName, const objectRegistry&);

        //- Remove named surface from specified registry
        static bool Delete(const word& surfName, const objectRegistry&);


    // Member Functions

    // Resolve iterator ambiguity in favour of Patch (not registry)

        using MeshReference::end;
        using MeshReference::cend;
        using MeshReference::begin;
        using MeshReference::cbegin;


    // Access

        //- Return the number of points
        virtual label nPoints() const
        {
            return MeshReference::points().size();
        }

        //- Return the number of faces
        virtual label nFaces() const
        {
            return MeshReference::size();
        }

        //- Return number of faces
        virtual label size() const
        {
            return MeshReference::size();
        }


        //- Return points
        virtual const pointField& points() const
        {
            return MeshReference::points();
        }

        //- Return faces
        virtual const faceList& faces() const
        {
            return this->storedFaces();
        }

        //- Const access to per-face zone/region information (demand-driven)
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
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


    // Modification

        //- Update with new contents.
        //- Removes existing fields if sizes have changed
        void copySurface
        (
            const pointField& points,
            const faceList& faces,
            //! validate the zone coverage (ignored)
            bool unused=false
        );

        //- Update with new contents.
        //- Removes existing fields if sizes have changed
        void copySurface
        (
            const meshedSurf& surf,
            //! validate the zone coverage (ignored)
            bool unused=false
        );

        //- Update with new contents.
        //- Removes existing fields if sizes have changed
        void copySurface
        (
            const MeshedSurface<face>& surf,
            //! validate the zone coverage (ignored)
            bool unused=false
        );

        //- Transfer the contents (and annul the parameters).
        //- Removes existing fields.
        void transfer
        (
            pointField&& points,
            faceList&& faces,
            labelList&& zoneIds = labelList()
        );

        //- Transfer the contents (and annul the parameters).
        //- Removes existing fields.
        void transfer
        (
            MeshedSurface<face>& surf,
            //! validate the zone coverage
            bool validate=false
        );


    // Fields

        //- Number of main entries, without PointData sub-registry
        label nFaceData() const;

        //- Number of entries on PointData sub-registry (if it exists)
        label nPointData() const;

        //- Query the field association (FACE or POINT)
        FieldAssociation queryFieldAssociation(const word& fieldName) const;

        //- Find the field object with the given name and required
        //- FieldAssociation (FACE or POINT).
        //  For FACE_POINT_DATA, face data are checked first.
        //  \return nullptr is the field was not found
        const regIOobject* findFieldObject
        (
            const word& fieldName,
            const FieldAssociation association
        ) const;

        //- General finding of the field object (FACE or POINT)
        //  Later specializations are used to restrict the scope.
        //  \return nullptr is the field was not found
        template<class GeoMeshType = void>
        const regIOobject* findFieldObject(const word& fieldName) const;

        //- General finding of the registry with the field object
        //- (FACE or POINT).
        //  Later specializations are used to restrict the scope.
        //  \return nullptr is the field was not found
        template<class GeoMeshType = void>
        const objectRegistry* whichRegistry(const word& fieldName) const;

        //- Regular data are stored directly on the registry
        template<class DataType>
        inline const objectRegistry& fieldData() const
        {
            return static_cast<const objectRegistry&>(*this);
        }

        //- Face data are stored directly on the registry
        const objectRegistry& faceData() const;

        //- Point data are stored in a sub-registry
        //  Note that this method will automatically create the corresponding
        //  sub-registry if it did not previously exist.
        //  Use the nPointData() methods instead if you wish to test for
        //  content without this side-effect.
        const objectRegistry& pointData() const;


        //- Create/store named zero field as face or point data
        //- (template parameter).
        //
        //  - Default is face-data (polySurfaceGeoMesh as template).
        //  - For point-data use polySurfacePointGeoMesh as template.
        //
        //  \return reference to the field
        template<class Type, class GeoMeshType = polySurfaceGeoMesh>
        DimensionedField<Type, GeoMeshType>&
        newField
        (
            const word& fieldName,
            const dimensionSet& dims
        );

        //- Copy/store named field as face or point data (template parameter).
        //
        //  - Default is face-data (polySurfaceGeoMesh as template).
        //  - For point-data use polySurfacePointGeoMesh as template.
        //
        //  \return reference to the field
        template<class Type, class GeoMeshType = polySurfaceGeoMesh>
        DimensionedField<Type, GeoMeshType>&
        storeField
        (
            const word& fieldName,
            const dimensionSet& dims,
            const Field<Type>& values
        );

        //- Move/store named field as face or point data (template parameter).
        //
        //  - Default is face-data (polySurfaceGeoMesh as template).
        //  - For point-data use polySurfacePointGeoMesh as template.
        //
        //  \return reference to the field
        template<class Type, class GeoMeshType = polySurfaceGeoMesh>
        DimensionedField<Type, GeoMeshType>&
        storeField
        (
            const word& fieldName,
            const dimensionSet& dims,
            Field<Type>&& values
        );


    // Writing

        //- Write - this is a no-op
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;


    // Storage Management

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
};


//- Find face field object (on main registry).
template<>
const regIOobject* polySurface::findFieldObject<polySurfaceGeoMesh>
(
    const word& fieldName
) const;


//- Find point field object (on sub-registry)
template<>
const regIOobject* polySurface::findFieldObject<polySurfacePointGeoMesh>
(
    const word& fieldName
) const;


//- Return field object (on main registry).
//- Face data field found on main registry?
template<>
const objectRegistry* polySurface::whichRegistry<polySurfaceGeoMesh>
(
    const word& fieldName
) const;


//- Point data field found on sub-registry?
template<>
const objectRegistry* polySurface::whichRegistry<polySurfaceGeoMesh>
(
    const word& fieldName
) const;


//- Point data are stored in a sub-registry
template<>
inline const objectRegistry&
polySurface::fieldData<polySurfacePointGeoMesh>() const
{
    return pointData();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "polySurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
