/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::sampledSurface

Group
    grpUtilitiesFunctionObjects

Description
    An abstract class for surfaces with sampling.

    The constructors for the derived classes should generally start in a
    'expired' condition (ie, needsUpdate() == true) and rely on a
    subsequent call to the update() method to complete the initialization.
    Delaying the final construction as late as possible allows the
    construction of surfaces that may depend on intermediate calculation
    results (eg, iso-surfaces) and also avoids the unnecessary
    reconstruction of surfaces between sampling intervals.

    It is the responsibility of the caller to ensure that the surface
    update() is called before the surface is used.  The update() method
    implementation should do nothing when the surface is already
    up-to-date.

    Any sampler is assumed to work for the standard volume field types.
    Some may also support surface fields.

    Dictionary entries:
    \table
        Property    | Description                           | Required | Default
        name        | Alternative name                      | no  |
        enabled     | Enable/disable the surface?           | no  | yes
        interpolate | Interpolate to nodes instead of faces | no  | false
        invariant   | Invariant with geometry change (use with caution!) | no  | false
    \endtable

Note
    The invariant switch is an advanced feature to declare that the surface
    is unaffected by changes in the general mesh geometry. For example, if
    sampling on a static patch while some other motion occurs elsewhere. If
    used improperly, there is a significant possibility for problems
    (caveat emptor).

SourceFiles
    sampledSurface.C
    sampledSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledSurface_H
#define Foam_sampledSurface_H

#include "polySurface.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMesh.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledSurface Declaration
\*---------------------------------------------------------------------------*/

class sampledSurface
:
    public meshedSurf
{
    // Private Data

        //- The name of the sample surface
        word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Should surface sampling be enabled?
        bool enabled_;

        //- Geometry is invariant (never changes)
        bool invariant_;

        //- Is point vs cell data
        bool isPointData_;

        //- Total surface area (demand-driven)
        mutable scalar area_;


protected:

    // Protected Member Functions

        //- Loop for sampling volume elements to faces.
        //  The defaultValue is used for invalid (negative) elements
        template<class Type>
        static tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler,
            const labelUList& elements,
            const faceList& fcs,
            const pointField& pts,
            const Type& defaultValue = Type(Zero)
        );


        //- Loop for interpolating volume elements to face points.
        template<class Type>
        static tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator,
            const labelUList& elements,
            const faceList& fcs,
            const pointField& pts
        );


        //- Create cell values by averaging the point values
        template<class Type>
        static tmp<VolumeField<Type>> pointAverage
        (
            const PointField<Type>& pfld
        );


        //- Additional cleanup when clearing the geometry
        virtual void clearGeom() const;

        //- Construct null
        explicit sampledSurface(const word& name, std::nullptr_t);


public:

    //- Runtime type information
    TypeName("sampledSurface");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        sampledSurface,
        word,
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        ),
        (name, mesh, dict)
    );


        //- PtrList read-construction helper
        class iNew
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<sampledSurface> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return sampledSurface::New(name, mesh_, dict);
            }
        };


        //- PtrList read-construction helper that captures dictionaries used
        //- during creation.
        class iNewCapture
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;

            //- Captured (recorded) dictionaries
            DynamicList<dictionary>& capture_;

        public:

            iNewCapture(const polyMesh& mesh, DynamicList<dictionary>& capture)
            :
                mesh_(mesh),
                capture_(capture)
            {}

            autoPtr<sampledSurface> operator()(Istream& is) const
            {
                word name(is);
                dictionary& dict = capture_.emplace_back(is);

                return sampledSurface::New(name, mesh_, dict);
            }
        };


    // Constructors

        //- Construct from name, mesh
        sampledSurface
        (
            const word& name,
            const polyMesh& mesh,
            const bool interpolateToPoints = false
        );

        //- Construct from dictionary
        sampledSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Clone
        autoPtr<sampledSurface> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a reference to the selected surface
        static autoPtr<sampledSurface> New
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor - calls clearGeom()
    virtual ~sampledSurface();


    // Member Functions

    // Access

        //- Access to the underlying mesh
        const polyMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Name of surface
        const word& name() const noexcept
        {
            return name_;
        }

        //- Surface is enabled
        bool enabled() const noexcept
        {
            return enabled_;
        }

        //- Surface is invariant with geometry change (caution)
        bool invariant() const noexcept
        {
            return invariant_;
        }

        //- Using interpolation to surface points
        bool isPointData() const noexcept
        {
            return isPointData_;
        }

        //- Change point/cell representation, may trigger an expire().
        //  \return old value
        virtual bool isPointData(const bool on);

        //- Does the surface need an update?
        virtual bool needsUpdate() const = 0;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire() = 0;

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was required
        virtual bool update() = 0;

        //- Points of surface
        virtual const pointField& points() const = 0;

        //- Faces of surface
        virtual const faceList& faces() const = 0;

        //- Face area vectors
        virtual const vectorField& Sf() const = 0;

        //- Face area magnitudes
        virtual const scalarField& magSf() const = 0;

        //- Face centres
        virtual const vectorField& Cf() const = 0;

        //- The total surface area
        scalar area() const;

        //- If element ids/order of the original surface are available
        virtual bool hasFaceIds() const
        {
            return false;
        }


    // General registry storage (optional)

        //- Get surface from registry if available.
        //  \param obr The objectRegistry to use
        //  \param lookupName Optional lookup name, use surface name if empty
        //  \return surface or nullptr
        polySurface* getRegistrySurface
        (
            const objectRegistry& obr,
            word lookupName = ""
        ) const;

        //- Copy surface into registry.
        //  \param obr The objectRegistry to use
        //  \param lookupName Optional lookup name, use surface name if empty
        //  \return surface or nullptr it surface should not be stored
        polySurface* storeRegistrySurface
        (
            objectRegistry& obr,
            word lookupName = ""
        ) const;

        //- Remove surface from registry.
        //  \param obr The objectRegistry to use
        //  \param lookupName Optional lookup name, use surface name if empty
        //  \return True if surface existed and was removed
        bool removeRegistrySurface
        (
            objectRegistry& obr,
            word lookupName = ""
        ) const;

        //- Copy/store sampled field onto registered surface (if it exists)
        template<class Type, class GeoMeshType>
        bool storeRegistryField
        (
            const objectRegistry& obr,
            const word& fieldName,
            const dimensionSet& dims,
            const Field<Type>& values,
            word lookupName = ""
        ) const;

        //- Move/store sampled field onto registered surface (if it exists)
        template<class Type, class GeoMeshType>
        bool storeRegistryField
        (
            const objectRegistry& obr,
            const word& fieldName,
            const dimensionSet& dims,
            Field<Type>&& values,
            word lookupName = ""
        ) const;


    // Sample (faces)

        //- Sample volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const = 0;

        //- Sample volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const = 0;

        //- Sample volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const = 0;

        //- Sample volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const = 0;

        //- Sample volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const = 0;


        //- Can it sample surface-fields?
        virtual bool withSurfaceFields() const;


        //- Sample surface field onto surface
        virtual tmp<scalarField> sample
        (
            const surfaceScalarField& sField
        ) const;

        //- Sample surface field onto surface
        virtual tmp<vectorField> sample
        (
            const surfaceVectorField& sField
        ) const;

        //- Sample surface field onto surface
        virtual tmp<sphericalTensorField> sample
        (
            const surfaceSphericalTensorField& sField
        ) const;

        //- Sample surface field onto surface
        virtual tmp<symmTensorField> sample
        (
            const surfaceSymmTensorField& sField
        ) const;

        //- Sample surface field onto surface
        virtual tmp<tensorField> sample
        (
            const surfaceTensorField& sField
        ) const;


    // Interpolate (points)

        //- Interpolate volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const = 0;

        //- Interpolate volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const = 0;

        //- Interpolate volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const = 0;

        //- Interpolate volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const = 0;

        //- Interpolate volume field onto surface points
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>& interpolator
        ) const = 0;


    // Edit

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }


    // Write

        //- Print information
        virtual void print(Ostream& os, int level=0) const;


    // Housekeeping

        //- Same as isPointData()
        bool interpolate() const noexcept { return isPointData_; }
};


// Global Operators

//- Ostream operator
Ostream& operator<<(Ostream& os, const sampledSurface& s);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
