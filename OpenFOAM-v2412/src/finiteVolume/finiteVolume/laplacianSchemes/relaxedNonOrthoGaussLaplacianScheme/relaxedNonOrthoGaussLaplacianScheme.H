/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::fv::relaxedNonOrthoGaussLaplacianScheme

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    laplacianSchemes
    {
        laplacian(<term>)    relaxedNonOrthoGauss <other options>;
    }
    \endverbatim

    and by using \c system/fvSolution:
    \verbatim
    relaxationFactors
    {
        equations
        {
            <term>           <relaxation factor>;
        }
    }
    \endverbatim

SourceFiles
    relaxedNonOrthoGaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef relaxedNonOrthoGaussLaplacianScheme_H
#define relaxedNonOrthoGaussLaplacianScheme_H

#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
             Class relaxedNonOrthoGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class relaxedNonOrthoGaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- No copy construct
        relaxedNonOrthoGaussLaplacianScheme
        (
            const relaxedNonOrthoGaussLaplacianScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const relaxedNonOrthoGaussLaplacianScheme&) = delete;


public:

    //- Runtime type information
    TypeName("relaxedNonOrthoGauss");


    // Constructors

        //- Construct null
        relaxedNonOrthoGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        relaxedNonOrthoGaussLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        relaxedNonOrthoGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    //- Destructor
    virtual ~relaxedNonOrthoGaussLaplacianScheme() = default;


    // Member Functions

        static tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& deltaCoeffs,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


// Use macros to emulate partial-specialisation of the Laplacian functions
// for scalar diffusivity gamma

#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
tmp<fvMatrix<Type>>                                                            \
relaxedNonOrthoGaussLaplacianScheme<Type, scalar>::fvmLaplacian                \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<GeometricField<Type, fvPatchField, volMesh>>                               \
relaxedNonOrthoGaussLaplacianScheme<Type, scalar>::fvcLaplacian                \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);


defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "relaxedNonOrthoGaussLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
