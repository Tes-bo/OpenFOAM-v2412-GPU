/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::sampledFaceZone

Description
    A sampledSurface defined by a faceZone or faceZones.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    faceZones;
            zones   (zone1 "sides.*");
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | faceZones                               | yes |
        zones    | zone selection as word/regex list       | yes |
        triangulate | triangulate faces                    | no  | false
    \endtable

SourceFiles
    sampledFaceZone.C
    sampledFaceZoneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledFaceZone_H
#define Foam_sampledFaceZone_H

#include "sampledSurface.H"
#include "MeshedSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledFaceZone Declaration
\*---------------------------------------------------------------------------*/

class sampledFaceZone
:
    public meshedSurface,
    public sampledSurface
{
    //- Mesh storage type
    typedef meshedSurface Mesh;


    // Private Data

        //- Selection (word/regex) of face zones
        const wordRes selectionNames_;

        //- IDs for selected face zones (sorted)
        mutable labelList zoneIds_;

        //- Triangulated faces or keep faces as is
        bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face. Is -1 for internal face
        labelList facePatchId_;


    // Private Member Functions

        //- Sample volume/boundary field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Sample surface field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const SurfaceField<Type>& sField
        ) const;

        //- Interpolate volume/boundary field onto surface points
        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;


        //- Re-map action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


protected:

        //- The selected face zones (sorted)
        const labelList& zoneIDs() const;


public:

    //- Runtime type information
    TypeName("faceZone");


    // Constructors

        //- Construct from components
        sampledFaceZone
        (
            const word& name,
            const polyMesh& mesh,
            const UList<wordRe>& zoneNames,
            const bool triangulate = false
        );

        //- Construct from dictionary
        sampledFaceZone
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledFaceZone() = default;


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Points of surface
        virtual const pointField& points() const
        {
            return Mesh::points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return Mesh::surfFaces();
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return labelList::null();
        }

        //- Face area vectors (normals)
        virtual const vectorField& Sf() const
        {
            return Mesh::Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return Mesh::magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return Mesh::Cf();
        }


    // Sample

        //- Sample volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const;

        //- Sample volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const;


        //- Can it sample surface-fields?
        virtual bool withSurfaceFields() const;


        //- Sample surface field on face zone
        virtual tmp<scalarField> sample
        (
            const surfaceScalarField&
        ) const;

        //- Sample surface field on face zone
        virtual tmp<vectorField> sample
        (
            const surfaceVectorField&
        ) const;

        //- Sample surface field on face zone
        virtual tmp<sphericalTensorField> sample
        (
            const surfaceSphericalTensorField&
        ) const;

        //- Sample surface field on face zone
        virtual tmp<symmTensorField> sample
        (
            const surfaceSymmTensorField&
        ) const;

        //- Sample surface field on face zone
        virtual tmp<tensorField> sample
        (
            const surfaceTensorField&
        ) const;


    // Interpolate

        //- Interpolate volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const;

        //- Interpolate volume field onto surface points
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>& interpolator
        ) const;


    // Output

        //- Print information
        virtual void print(Ostream& os, int level=0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledFaceZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
