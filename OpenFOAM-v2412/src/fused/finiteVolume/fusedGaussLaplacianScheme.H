/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2024 M. Janssens
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
    Foam::fv::GaussLaplacianScheme

Group
    grpFvLaplacianSchemes

Description
    Variant of GaussLaplacian that avoids intermediate fields

SourceFiles
    fusedGaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef fusedGaussLaplacianScheme_H
#define fusedGaussLaplacianScheme_H

#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class fusedGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class fusedGaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        //- See gaussLaplacianScheme
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- No copy construct
        fusedGaussLaplacianScheme(const fusedGaussLaplacianScheme&)
            = delete;

        //- No copy assignment
        void operator=(const fusedGaussLaplacianScheme&) = delete;


public:

    //- Runtime type information
    TypeName("fusedGauss");


    // Constructors

        //- Construct null
        fusedGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        fusedGaussLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        fusedGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    //- Destructor
    virtual ~fusedGaussLaplacianScheme() = default;


    // Member Functions

        static tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        ////- Helper: grad on component. Move to GaussGrad? Template on
        ////  - cell-access operator?
        ////  - cell-to-face operator?
        ////  - face-to-cell operator?
        //static void gradComponent
        //(
        //    const surfaceScalarField& weights,
        //    const GeometricField<Type, fvPatchField, volMesh>& vf,
        //    const direction cmpt,
        //    GeometricField<Type, fvPatchField, volMesh>& gGrad
        //);

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // New: override laplacianScheme::fvc/mLaplacian for scalar gamma and
    //      scalar,vector volField

        virtual tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvPatchField, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// Use macros to emulate partial-specialisation of the Laplacian functions
// for scalar diffusivity gamma

#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
tmp<fvMatrix<Type>> fusedGaussLaplacianScheme<Type, scalar>::fvmLaplacian      \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<GeometricField<Type, fvPatchField, volMesh>>                               \
fusedGaussLaplacianScheme<Type, scalar>::fvcLaplacian                          \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);

defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);

// Unweighted laplacian
template<>
tmp<GeometricField<scalar, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<scalar, scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&
);
template<>
tmp<GeometricField<vector, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<vector, scalar>::fvcLaplacian
(
    const GeometricField<vector, fvPatchField, volMesh>&
);

// Weighted laplacian
template<>
tmp<GeometricField<scalar, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<scalar, scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&,
    const GeometricField<scalar, fvPatchField, volMesh>&
);
template<>
tmp<GeometricField<vector, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<vector, scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&,
    const GeometricField<vector, fvPatchField, volMesh>&
);
template<>
tmp<fvMatrix<scalar>>
fusedGaussLaplacianScheme<scalar, scalar>::fvmLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&,
    const GeometricField<scalar, fvPatchField, volMesh>&
);
template<>
tmp<fvMatrix<vector>>
fusedGaussLaplacianScheme<vector, scalar>::fvmLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>&,
    const GeometricField<vector, fvPatchField, volMesh>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fusedGaussLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
