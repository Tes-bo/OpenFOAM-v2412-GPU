/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::sampledPatchInternalField

Description
    Variation of sampledPatch that samples the internalField (at a given
    normal distance from the patch) instead of the patchField.
    Note:
    - interpolate=false : get cell value on faces
    - interpolate=true  : interpolate inside cell and interpolate to points
    There is no option to get interpolated value inside the cell on the faces.

    This is often embedded as part of a sampled surfaces function object.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type    patchInternalField;
            patches (inlet "outlet.*");
            offsetMode  normal;
            distance    0.05;
        }
    }
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | patchInternalField                      | yes      |
        patches  | patch selection as word/regex list      | yes      |
        offsetMode | normal/uniform/nonuniform             | no       | normal
        distance | distance for normal offset              | partly   |
        offset   | point offset for uniform offset         | partly   |
        offsets  | point offsets for nonuniform offset     | partly   |
    \endtable

SourceFiles
    sampledPatchInternalField.C

\*---------------------------------------------------------------------------*/

#ifndef sampledPatchInternalField_H
#define sampledPatchInternalField_H

#include "sampledPatch.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class sampledPatchInternalField Declaration
\*---------------------------------------------------------------------------*/

class sampledPatchInternalField
:
    public sampledPatch
{
    // Private data

        //- Mapping engines
        PtrList<mappedPatchBase> mappers_;


    // Private Member Functions

        //- Sample volume field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;


public:

    //- Runtime type information
    TypeName("sampledPatchInternalField");


    // Constructors

        //- Construct from dictionary
        sampledPatchInternalField
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledPatchInternalField() = default;


    // Member Functions

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
    #include "sampledPatchInternalFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
