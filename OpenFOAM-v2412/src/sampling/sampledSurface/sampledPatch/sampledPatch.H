/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::sampledPatch

Description
    A sampledSurface on patches. Non-triangulated by default.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    patch;
            patches (inlet "outlet.*");
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | patch                                   | yes      |
        patches  | patch selection as word/regex list      | yes      |
        triangulate | triangulate faces                    | no       | false
    \endtable

SourceFiles
    sampledPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledPatch_H
#define Foam_sampledPatch_H

#include "sampledSurface.H"
#include "MeshedSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledPatch Declaration
\*---------------------------------------------------------------------------*/

class sampledPatch
:
    public meshedSurface,
    public sampledSurface
{
    //- Mesh storage type
    typedef meshedSurface Mesh;


    // Private Data

        //- Selection (word/regex) of patches
        const wordRes selectionNames_;

        //- The IDs for selected patches. Sorted and no emptyPolyPatch
        mutable labelList patchIDs_;

        //- Triangulated faces or keep faces as is
        bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- Start indices (in patchFaceLabels_) of patches
        labelList patchStart_;

        //- For every face (or triangle) the originating patch
        labelList patchIndex_;

        //- For every face (or triangle) the index in the originating patch
        labelList patchFaceLabels_;


    // Private Member Functions

        //- Sample boundary field (from volume field) onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Sample boundary field (from surface field) onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const SurfaceField<Type>& sField
        ) const;

        //- Interpolate boundary field (from volume field) onto surface points
        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;


        //- Re-map action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


protected:

        //- The selection (word/regex) of patches
        const wordRes& patchNames() const
        {
            return selectionNames_;
        }

        //- The patches selected
        const labelList& patchIDs() const;

        //- The offset into patchIndex, patchFaceLabels
        const labelList& patchStart() const
        {
            return patchStart_;
        }

        //- For each face, the patch ID
        const labelList& patchIndex() const
        {
            return patchIndex_;
        }

        //- For each face, the patch local face ID
        const labelList& patchFaceLabels() const
        {
            return patchFaceLabels_;
        }


public:

    //- Runtime type information
    TypeName("sampledPatch");


    // Constructors

        //- Construct from components
        sampledPatch
        (
            const word& name,
            const polyMesh& mesh,
            const UList<wordRe>& patchNames,
            const bool triangulate = false
        );

        //- Construct from dictionary
        sampledPatch
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledPatch() = default;


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

        //- Face area vectors
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

        //- Sample boundary of volume field onto surface faces
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>& sampler
        ) const;

        //- Sample boundary of volume field onto surface faces
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>& sampler
        ) const;

        //- Sample boundary of volume field onto surface faces
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>& sampler
        ) const;

        //- Sample boundary of volume field onto surface faces
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>& sampler
        ) const;

        //- Sample boundary of volume field onto surface faces
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>& sampler
        ) const;


        //- Can it sample surface-fields?
        virtual bool withSurfaceFields() const;


        //- Sample boundary of surface field on surface
        virtual tmp<scalarField> sample
        (
            const surfaceScalarField&
        ) const;

        //- Sample boundary of surface field on surface
        virtual tmp<vectorField> sample
        (
            const surfaceVectorField&
        ) const;

        //- Sample boundary of surface field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const surfaceSphericalTensorField&
        ) const;

        //- Sample boundary of surface field on surface
        virtual tmp<symmTensorField> sample
        (
            const surfaceSymmTensorField&
        ) const;

        //- Sample boundary of surface field on surface
        virtual tmp<tensorField> sample
        (
            const surfaceTensorField&
        ) const;


    // Interpolate

        //- Interpolate boundary of volume field onto surface points
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>& interpolator
        ) const;

        //- Interpolate boundary of volume field onto surface points
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>& interpolator
        ) const;

        //- Interpolate boundary of volume field onto surface points
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>& interpolator
        ) const;

        //- Interpolate boundary of volume field onto surface points
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>& interpolator
        ) const;

        //- Interpolate boundary of volume field onto surface points
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
    #include "sampledPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
