/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fusedGaussLaplacianScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvLaplacianScheme(fusedGaussLaplacianScheme)

#define declareFvmLaplacianScalarGamma(Type)                                   \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvMatrix<Foam::Type>>                                          \
Foam::fv::fusedGaussLaplacianScheme<Foam::Type, Foam::scalar>::                \
fvmLaplacian                                                                   \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    DebugPout<< "fusedGaussLaplacianScheme::fvmLaplacian on " << vf.name()     \
        << " with scalar gamma " << gamma.name() << endl;                      \
                                                                               \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf              \
    (                                                                          \
        gamma*mesh.magSf()                                                     \
    );                                                                         \
                                                                               \
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected                         \
    (                                                                          \
        gammaMagSf,                                                            \
        this->tsnGradScheme_().deltaCoeffs(vf),                                \
        vf                                                                     \
    );                                                                         \
    fvMatrix<Type>& fvm = tfvm.ref();                                          \
                                                                               \
    if (this->tsnGradScheme_().corrected())                                    \
    {                                                                          \
        if (mesh.fluxRequired(vf.name()))                                      \
        {                                                                      \
            fvm.faceFluxCorrectionPtr() = std::make_unique                     \
            <                                                                  \
                GeometricField<Type, fvsPatchField, surfaceMesh>               \
            >                                                                  \
            (                                                                  \
                gammaMagSf*this->tsnGradScheme_().correction(vf)               \
            );                                                                 \
                                                                               \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    *fvm.faceFluxCorrectionPtr()                               \
                )().primitiveField();                                          \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)           \
                )().primitiveField();                                          \
        }                                                                      \
    }                                                                          \
                                                                               \
    return tfvm;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::GeometricField<Foam::Type, Foam::fvPatchField, Foam::volMesh>> \
Foam::fv::fusedGaussLaplacianScheme<Foam::Type, Foam::scalar>::                \
fvcLaplacian                                                                   \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& gamma,           \
    const GeometricField<Type, fvPatchField, volMesh>& vf                      \
)                                                                              \
{                                                                              \
    DebugPout<< "fvcLaplacian on " << vf.name()                                \
        << " with scalar gamma " << gamma.name() << endl;                      \
                                                                               \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian                \
    (                                                                          \
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())         \
    );                                                                         \
                                                                               \
    tLaplacian.ref().rename                                                    \
    (                                                                          \
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'                    \
    );                                                                         \
                                                                               \
    return tLaplacian;                                                         \
}


template<>
Foam::tmp
<
    Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>
>
Foam::fv::fusedGaussLaplacianScheme<Foam::scalar, Foam::scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>& gamma,
    const GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    typedef scalar Type;
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    DebugPout
        << "fusedGaussLaplacianScheme<scalar, scalar>::fvcLaplacian"
        << " on " << vf.name() << " with gamma " << gamma.name() << endl;

    const fvMesh& mesh = vf.mesh();

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                gamma.dimensions()*vf.dimensions()/dimArea, Zero
            ),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();

    const auto tweights(this->tinterpGammaScheme_().weights(gamma));
    const auto& weights = tweights();
    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& deltaCoeffs = tdeltaCoeffs();

    if (this->tsnGradScheme_().corrected())
    {
        // Calculate sn gradient
        tmp<SurfaceFieldType> tfaceGrad
        (
            new SurfaceFieldType
            (
                IOobject
                (
                    "snGradCorr("+vf.name()+')',
                    vf.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                vf.dimensions()
            )
        );

        {
            // Calculate gradient
            tmp<GradFieldType> tgGrad
            (
                gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme("grad(" + vf.name() + ')')
                )().grad(vf, "grad(" + vf.name() + ')')
            );
            const auto& gGrad = tgGrad();

            // Doing a dotinterpolate with nonOrthCorrectionVectors
            const auto dotInterpolate = [&]
            (
                const vector& area,
                const scalar lambda,

                const GradType& ownVal,
                const GradType& neiVal,

                const vector& dotVector,

                Type& result
            )
            {
                result = dotVector&(lambda*(ownVal - neiVal) + neiVal);
            };

            fvc::interpolate
            (
                mesh.surfaceInterpolation::weights(),   // linear interpolation
                gGrad,                          // volume field
                mesh.nonOrthCorrectionVectors(),// surface multiplier
                dotInterpolate,
                tfaceGrad.ref()
            );
        }
        const auto& faceGrad = tfaceGrad();

        const auto snGrad = [&]
        (
            const vector& Sf,

            const scalar weight,
            const scalar ownGamma,
            const scalar neiGamma,

            const scalar dc,
            const Type& ownVal,
            const Type& neiVal,

            const Type& correction
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal) + correction);
            const scalar faceGamma(weight*(ownGamma-neiGamma)+neiGamma);
            return mag(Sf)*faceGamma*snGrad;
        };

        fvc::surfaceSnSum
        (
            weights,    // gamma weights
            gamma,

            deltaCoeffs,
            vf,

            faceGrad,   // face-based addition

            snGrad,

            result,
            false       // avoid boundary evaluation until volume division
        );
    }
    else
    {
        const auto snGrad = [&]
        (
            const vector& Sf,

            const scalar weight,
            const scalar ownGamma,
            const scalar neiGamma,

            const scalar dc,
            const Type& ownVal,
            const Type& neiVal
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal));
            const scalar faceGamma(weight*(ownGamma-neiGamma)+neiGamma);
            return mag(Sf)*faceGamma*snGrad;
        };

        fvc::surfaceSnSum
        (
            weights,
            gamma,

            deltaCoeffs,
            vf,

            snGrad,

            result,
            false       // avoid boundary evaluation until volume division
        );
    }

    result.primitiveFieldRef() /= mesh.V();
    result.correctBoundaryConditions();

    return tresult;
}


template<>
Foam::tmp
<
    Foam::GeometricField<Foam::vector, Foam::fvPatchField, Foam::volMesh>
>
Foam::fv::fusedGaussLaplacianScheme<Foam::vector, Foam::scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>& gamma,
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    DebugPout
        << "fusedGaussLaplacianScheme<vector, scalar>::fvcLaplacian"
        << " on " << vf.name() << " with gamma " << gamma.name() << endl;

    typedef vector Type;
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = vf.mesh();

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>
            (
                gamma.dimensions()*vf.dimensions()/dimArea, Zero
            ),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();

    const auto tweights(this->tinterpGammaScheme_().weights(gamma));
    const auto& weights = tweights();
    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& deltaCoeffs = tdeltaCoeffs();

    if (this->tsnGradScheme_().corrected())
    {
        // Calculate sn gradient
        tmp<SurfaceFieldType> tfaceGrad
        (
            new SurfaceFieldType
            (
                IOobject
                (
                    "snGradCorr("+vf.name()+')',
                    vf.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                vf.dimensions()
            )
        );

        {
            // Calculate gradient
            tmp<GradFieldType> tgGrad
            (
                gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme("grad(" + vf.name() + ')')
                )().grad(vf, "grad(" + vf.name() + ')')
            );
            const auto& gGrad = tgGrad();

            // Doing a dotinterpolate with nonOrthCorrectionVectors
            const auto dotInterpolate = [&]
            (
                const vector& area,
                const scalar lambda,

                const GradType& ownVal,
                const GradType& neiVal,

                const vector& dotVector,

                Type& result
            )
            {
                result = dotVector&(lambda*(ownVal - neiVal) + neiVal);
            };

            fvc::interpolate
            (
                mesh.surfaceInterpolation::weights(),   // linear interpolation
                gGrad,                          // volume field
                mesh.nonOrthCorrectionVectors(),// surface multiplier
                dotInterpolate,
                tfaceGrad.ref()
            );
        }
        const auto& faceGrad = tfaceGrad();

        const auto snGrad = [&]
        (
            const vector& Sf,

            const scalar weight,
            const scalar ownGamma,
            const scalar neiGamma,

            const scalar dc,
            const Type& ownVal,
            const Type& neiVal,

            const Type& correction
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal) + correction);
            const scalar faceGamma(weight*(ownGamma-neiGamma)+neiGamma);
            return mag(Sf)*faceGamma*snGrad;
        };

        fvc::surfaceSnSum
        (
            weights,    // gamma weights
            gamma,

            deltaCoeffs,
            vf,

            faceGrad,   // face-based addition

            snGrad,

            result,
            false       // avoid boundary evaluation until volume division
        );
    }
    else
    {
        const auto snGrad = [&]
        (
            const vector& Sf,

            const scalar weight,
            const scalar ownGamma,
            const scalar neiGamma,

            const scalar dc,
            const Type& ownVal,
            const Type& neiVal
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal));
            const scalar faceGamma(weight*(ownGamma-neiGamma)+neiGamma);
            return mag(Sf)*faceGamma*snGrad;
        };

        fvc::surfaceSnSum
        (
            weights,    // gamma weights
            gamma,

            deltaCoeffs,
            vf,

            snGrad,

            result,
            false       // avoid boundary evaluation until volume division
        );
    }

    result.primitiveFieldRef() /= mesh.V();
    result.correctBoundaryConditions();

    return tresult;
}


template<>
Foam::tmp<Foam::fvMatrix<Foam::scalar>>
Foam::fv::fusedGaussLaplacianScheme<Foam::scalar, Foam::scalar>::fvmLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>& gamma,
    const GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    // TBD
    DebugPout
        << "fusedGaussLaplacianScheme<scalar, scalar>::fvmLaplacian"
        << " on " << vf.name() << " with gamma " << gamma.name() << endl;
    return fvmLaplacian(this->tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<>
Foam::tmp<Foam::fvMatrix<Foam::vector>>
Foam::fv::fusedGaussLaplacianScheme<Foam::vector, Foam::scalar>::fvmLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>& gamma,
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    // TBD
    DebugPout
        << "fusedGaussLaplacianScheme<vector, scalar>::fvmLaplacian"
        << " on " << vf.name() << " with gamma " << gamma.name() << endl;
    return fvmLaplacian(this->tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<>
Foam::tmp
<
    Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>
>
Foam::fv::fusedGaussLaplacianScheme<Foam::scalar, Foam::scalar>::fvcLaplacian
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
    typedef scalar Type;

    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = vf.mesh();

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(vf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();

    DebugPout
        << "fusedGaussLaplacianScheme<scalar, GType>::fvcLaplacian on "
        << vf.name()
        << " to generate " << result.name() << endl;


    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& deltaCoeffs = tdeltaCoeffs();


    if (this->tsnGradScheme_().corrected())
    {
        // Calculate sn gradient
        tmp<SurfaceFieldType> tfaceGrad
        (
            new SurfaceFieldType
            (
                IOobject
                (
                    "snGradCorr("+vf.name()+')',
                    vf.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                vf.dimensions()
            )
        );

        {
            // Calculate gradient
            tmp<GradFieldType> tgGrad
            (
                gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme("grad(" + vf.name() + ')')
                )().grad(vf, "grad(" + vf.name() + ')')
            );
            const auto& gGrad = tgGrad();

            // Doing a dotinterpolate with nonOrthCorrectionVectors
            const auto dotInterpolate = [&]
            (
                const vector& area,
                const scalar lambda,

                const GradType& ownVal,
                const GradType& neiVal,

                const vector& dotVector,

                Type& result
            )
            {
                result = dotVector&(lambda*(ownVal - neiVal) + neiVal);
            };

            fvc::interpolate
            (
                mesh.surfaceInterpolation::weights(),   // linear interpolation
                gGrad,                          // volume field
                mesh.nonOrthCorrectionVectors(),// surface multiplier
                dotInterpolate,
                tfaceGrad.ref()
            );
        }
        const auto& faceGrad = tfaceGrad();

        const auto snGrad = [&]
        (
            const vector& Sf,
            const scalar dc,
            const Type& ownVal,
            const Type& neiVal,
            const Type& correction
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal) + correction);
            return mag(Sf)*snGrad;
        };

        fvc::surfaceSnSum
        (
            deltaCoeffs,
            vf,
            faceGrad,   // face-based addition
            snGrad,
            result,
            false       // avoid boundary evaluation until volume division
        );
    }
    else
    {
        const auto snGrad = [&]
        (
            const vector& Sf,
            const scalar dc,
            const Type& ownVal,
            const Type& neiVal
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal));
            return mag(Sf)*snGrad;
        };

        fvc::surfaceSnSum
        (
            deltaCoeffs,
            vf,
            snGrad,
            result,
            false       // avoid boundary evaluation until volume division
        );
    }

    result.primitiveFieldRef() /= mesh.V();
    result.correctBoundaryConditions();

    return tresult;
}


template<>
Foam::tmp
<
    Foam::GeometricField<Foam::vector, Foam::fvPatchField, Foam::volMesh>
>
Foam::fv::fusedGaussLaplacianScheme<Foam::vector, Foam::scalar>::fvcLaplacian
(
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    typedef vector Type;

    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = vf.mesh();

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(vf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();

    DebugPout
        << "fusedGaussLaplacianScheme<vector, GType>::fvcLaplacian on "
        << vf.name()
        << " to generate " << result.name() << endl;


    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& deltaCoeffs = tdeltaCoeffs();


    if (this->tsnGradScheme_().corrected())
    {
        // Calculate sn gradient
        tmp<SurfaceFieldType> tfaceGrad
        (
            new SurfaceFieldType
            (
                IOobject
                (
                    "snGradCorr("+vf.name()+')',
                    vf.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                vf.dimensions()
            )
        );

        {
            // Calculate gradient
            tmp<GradFieldType> tgGrad
            (
                gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme("grad(" + vf.name() + ')')
                )().grad(vf, "grad(" + vf.name() + ')')
            );
            const auto& gGrad = tgGrad();

            // Doing a dotinterpolate with nonOrthCorrectionVectors
            const auto dotInterpolate = [&]
            (
                const vector& area,
                const scalar lambda,

                const GradType& ownVal,
                const GradType& neiVal,

                const vector& dotVector,

                Type& result
            )
            {
                result = dotVector&(lambda*(ownVal - neiVal) + neiVal);
            };

            fvc::interpolate
            (
                mesh.surfaceInterpolation::weights(),   // linear interpolation
                gGrad,                          // volume field
                mesh.nonOrthCorrectionVectors(),// surface multiplier
                dotInterpolate,
                tfaceGrad.ref()
            );
        }
        const auto& faceGrad = tfaceGrad();

        const auto snGrad = [&]
        (
            const vector& Sf,
            const scalar dc,
            const Type& ownVal,
            const Type& neiVal,
            const Type& correction
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal) + correction);
            return mag(Sf)*snGrad;
        };

        fvc::surfaceSnSum
        (
            deltaCoeffs,
            vf,
            faceGrad,   // face-based addition
            snGrad,
            result,
            false       // avoid boundary evaluation until volume division
        );
    }
    else
    {
        const auto snGrad = [&]
        (
            const vector& Sf,
            const scalar dc,
            const Type& ownVal,
            const Type& neiVal
        ) -> Type
        {
            const auto snGrad(dc*(neiVal-ownVal));
            return mag(Sf)*snGrad;
        };

        fvc::surfaceSnSum
        (
            deltaCoeffs,
            vf,
            snGrad,
            result,
            false       // avoid boundary evaluation until volume division
        );
    }

    result.primitiveFieldRef() /= mesh.V();
    result.correctBoundaryConditions();

    return tresult;
}


declareFvmLaplacianScalarGamma(scalar);
declareFvmLaplacianScalarGamma(vector);
declareFvmLaplacianScalarGamma(sphericalTensor);
declareFvmLaplacianScalarGamma(symmTensor);
declareFvmLaplacianScalarGamma(tensor);


// ************************************************************************* //
