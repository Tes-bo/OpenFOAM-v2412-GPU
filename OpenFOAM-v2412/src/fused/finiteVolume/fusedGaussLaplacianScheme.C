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
#include "fvcSurfaceOps.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
fusedGaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DebugPout
        << "fusedGaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected on "
        << vf.name()
        << " with gammaMagSf " << gammaMagSf.name()
        << " with deltaCoeffs " << deltaCoeffs.name()
        << endl;

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    //fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
    multiply
    (
        fvm.upper(),
        deltaCoeffs.primitiveField(),
        gammaMagSf.primitiveField()
    );
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchi)
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
        const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];

        auto& intCoeffs = fvm.internalCoeffs()[patchi];
        auto& bouCoeffs = fvm.boundaryCoeffs()[patchi];

        if (pvf.coupled())
        {
            //intCoeffs = pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            multiply
            (
                intCoeffs,
                pGamma,
                pvf.gradientInternalCoeffs(pDeltaCoeffs)()
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            multiply
            (
                bouCoeffs,
                pGamma,
                pvf.gradientBoundaryCoeffs(pDeltaCoeffs)()
            );
            bouCoeffs.negate();
        }
        else
        {
            //intCoeffs = pGamma*pvf.gradientInternalCoeffs();
            multiply
            (
                intCoeffs,
                pGamma,
                pvf.gradientInternalCoeffs()()
            );
            //bouCoeffs = -pGamma*pvf.gradientBoundaryCoeffs();
            multiply
            (
                bouCoeffs,
                pGamma,
                pvf.gradientBoundaryCoeffs()()
            );
            bouCoeffs.negate();
        }
    }

    return tfvm;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
fusedGaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
    const surfaceVectorField& SfGammaCorr,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    DebugPout<< "fusedGaussLaplacianScheme<Type, GType>::gammaSnGradCorr on "
        << vf.name() << " with SfGammCorr " << SfGammaCorr.name() << endl;

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tgammaSnGradCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "gammaSnGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            SfGammaCorr.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    auto& gammaSnGradCorr = tgammaSnGradCorr.ref();
    gammaSnGradCorr.oriented() = SfGammaCorr.oriented();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        gammaSnGradCorr.replace
        (
            cmpt,
            fvc::dotInterpolate(SfGammaCorr, fvc::grad(vf.component(cmpt)))
        );
    }

    return tgammaSnGradCorr;
}


/*
template<class Type, class GType>
void fusedGaussLaplacianScheme<Type, GType>::gradComponent
(
    const surfaceScalarField& weights,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const direction cmpt,
    GeometricField<Type, fvPatchField, volMesh>& gGrad
)
{
    gGrad = dimensioned<Type>(vf.dimensions()/dimLength, Zero);

    // Calculate grad of vf.component(cmpt)
    fvc::GaussOp
    (
        vf,
        weights,
        componentInterpolate<Type>(cmpt),
        gGrad
    );
}
template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
fusedGaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
    const surfaceScalarField& weights,
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tgammaSnGradCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "gammaSnGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            gamma.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    auto& gammaSnGradCorr = tgammaSnGradCorr.ref();

    gammaSnGradCorr.oriented() = gamma.oriented();


    GeometricField<Type, fvPatchField, volMesh> gradCmptFld
    (
        IOobject
        (
            vf.name() + ".component()",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        vf.dimensions()
    );

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        // Calculate fvc::grad(vf.component(cmpt)) into gradCmptFld
        gradComponent
        (
            weights,
            vf,
            cmpt,

            gradCmptFld
        );

        //gammaSnGradCorr.replace
        //(
        //    cmpt,
        //    fvc::dotInterpolate(SfGammaCorr, gradCmptFld)
        //);


        fvc::interpolate
        (
            weights,

            gradCmptFld,        // fvc::grad(vf.component(cmpt))

            gamma,              // weight field

            tanInterpolate<Type, GType>(cmpt),

            gammaSnGradCorr
        );
    }

    return tgammaSnGradCorr;
}
*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    //typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    tmp<FieldType> tresult
    (
        new FieldType
        (
            IOobject
            (
                "laplacian(" + vf.name() + ')',
                vf.instance(),
                vf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            vf.mesh(),
            dimensioned<Type>(vf.dimensions()/dimArea, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    FieldType& result = tresult.ref();

    DebugPout
        << "fusedGaussLaplacianScheme<Type, GType>::fvcLaplacian on "
        << vf.name()
        << " to generate " << result.name() << endl;

    // Note: cannot use fvc::GaussOp since specialised handling on boundary.
    // Maybe bypass for processor boundaries?

    const auto tdeltaCoeffs(this->tsnGradScheme_().deltaCoeffs(vf));
    const auto& deltaCoeffs = tdeltaCoeffs();


    const fvMesh& mesh = vf.mesh();

    if (this->tsnGradScheme_().corrected())
    {
        // Problem when instantiating for tensors - outerProduct not defined.
        // scalar/vector specialisations in *Schemes.C file.

        FatalErrorInFunction<< "Corrected snGrad not supported for field "
            << vf.name() << exit(FatalError);

        //typedef typename outerProduct<vector, Type>::type GradType;
        //typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;
        //// Calculate sn gradient
        //tmp<SurfaceFieldType> tfaceGrad
        //(
        //    new SurfaceFieldType
        //    (
        //        IOobject
        //        (
        //            "snGradCorr("+vf.name()+')',
        //            vf.instance(),
        //            mesh,
        //            IOobject::NO_READ,
        //            IOobject::NO_WRITE
        //        ),
        //        mesh,
        //        vf.dimensions()
        //    )
        //);
        //
        //{
        //    // Calculate gradient
        //    tmp<GradFieldType> tgGrad
        //    (
        //        gradScheme<Type>::New
        //        (
        //            mesh,
        //            mesh.gradScheme("grad(" + vf.name() + ')')
        //        )().grad(vf, "grad(" + vf.name() + ')')
        //    );
        //    const auto& gGrad = tgGrad();
        //
        //    // Doing a dotinterpolate with nonOrthCorrectionVectors
        //    const auto dotInterpolate = [&]
        //    (
        //        const vector& area,
        //        const scalar lambda,
        //
        //        const GradType& ownVal,
        //        const GradType& neiVal,
        //
        //        const vector& dotVector,
        //
        //        Type& result
        //    )
        //    {
        //        result = dotVector&(lambda*(ownVal - neiVal) + neiVal);
        //    };
        //
        //    fvc::interpolate
        //    (
        //        mesh.surfaceInterpolation::weights(),   // linear interp
        //        gGrad,                          // volume field
        //        mesh.nonOrthCorrectionVectors(),// surface multiplier
        //        dotInterpolate,
        //        tfaceGrad.ref()
        //    );
        //}
        //const auto& faceGrad = tfaceGrad();
        //
        //const auto snGrad = [&]
        //(
        //    const vector& Sf,
        //    const scalar dc,
        //    const Type& ownVal,
        //    const Type& neiVal,
        //    const Type& correction
        //) -> Type
        //{
        //    const auto snGrad(dc*(neiVal-ownVal) + correction);
        //    return mag(Sf)*snGrad;
        //};
        //
        //fvc::surfaceSnSum
        //(
        //    deltaCoeffs,
        //    vf,
        //    faceGrad,   // face-based addition
        //    snGrad,
        //    result,
        //    false       // avoid boundary evaluation until volume division
        //);
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


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DebugPout
        << "fusedGaussLaplacianScheme<Type, GType>::fvcLaplacian on "
        << vf.name() << " with gamma " << gamma.name() << endl;

    return fvcLaplacian(this->tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
fusedGaussLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DebugPout<< "fusedGaussLaplacianScheme<Type, GType>::fvmLaplacian on "
        << vf.name() << " with gamma " << gamma.name() << endl;

    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected
    (
        SfGammaSn,
        this->tsnGradScheme_().deltaCoeffs(vf),
        vf
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfaceFluxCorrection
        = gammaSnGradCorr(SfGammaCorr, vf);

    if (this->tsnGradScheme_().corrected())
    {
        tfaceFluxCorrection.ref() +=
            SfGammaSn*this->tsnGradScheme_().correction(vf);
    }

    fvm.source() -= mesh.V()*fvc::div(tfaceFluxCorrection())().primitiveField();

    if (mesh.fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr(tfaceFluxCorrection.ptr());
    }

    return tfvm;
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
fusedGaussLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DebugPout<< "fusedGaussLaplacianScheme<Type, GType>::fvmLaplacian on "
        << vf.name() << " with gamma " << gamma.name() << endl;
    return fvmLaplacian(this->tinterpGammaScheme_().interpolate(gamma)(), vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
fusedGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DebugPout<< "fusedGaussLaplacianScheme<Type, GType>::fvcLaplacian on "
        << vf.name() << " with gamma " << gamma.name() << endl;

    const fvMesh& mesh = this->mesh();

    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);

    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
    (
        fvc::div
        (
            SfGammaSn*this->tsnGradScheme_().snGrad(vf)
          + gammaSnGradCorr(SfGammaCorr, vf)
        )
    );

    tLaplacian.ref().rename
    (
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
