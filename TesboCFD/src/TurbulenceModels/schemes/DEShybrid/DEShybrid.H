/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
    Copyright (C) 2022 Upstream CFD GmbH
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
    Foam::DEShybrid

Description
    Improved hybrid convection scheme of Travin et al. for hybrid RAS/LES
    calculations with enhanced Grey Area Mitigation (GAM) behaviour.

    The scheme provides a blend between two convection schemes, based on local
    properties including the wall distance, velocity gradient and eddy
    viscosity.  The scheme was originally developed for DES calculations to
    blend a low-dissipative scheme, e.g. linear, in the vorticity-dominated,
    finely-resolved regions and a numerically more robust, e.g. upwind-biased,
    convection scheme in irrotational or coarsely-resolved regions.

    The routine calculates the blending factor denoted as "sigma" in the
    literature reference, where 0 <= sigma <= sigmaMax, which is then employed
    to set the weights:
    \f[
        weight = (1-sigma) w_1 + sigma w_2
    \f]

    where
    \vartable
        sigma | blending factor
        w_1   | scheme 1 weights
        w_2   | scheme 2 weights
    \endvartable

    First published in:
    \verbatim
        Travin, A., Shur, M., Strelets, M., & Spalart, P. R. (2000).
        Physical and numerical upgrades in the detached-eddy
        simulation of complex turbulent flows.
        In LES of Complex Transitional and Turbulent Flows.
        Proceedings of the Euromech Colloquium 412. Munich, Germany
    \endverbatim

    Original publication contained a typo for \c C_H3 constant.
    Corrected version with minor changes for 2 lower limiters published in:
    \verbatim
        Spalart, P., Shur, M., Strelets, M., & Travin, A. (2012).
        Sensitivity of landing-gear noise predictions by large-eddy
        simulation to numerics and resolution.
        In 50th AIAA Aerospace Sciences Meeting Including the
        New Horizons Forum and Aerospace Exposition. Nashville, US.
        DOI:10.2514/6.2012-1174
    \endverbatim

    Example of the \c DEShybrid scheme specification using \c linear
    within the LES region and \c linearUpwind within the RAS region:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,U)      Gauss DEShybrid
            linear                    // scheme 1
            linearUpwind grad(U)      // scheme 2
            delta                     // LES delta name, e.g. 'delta', 'hmax'
            0.65                      // CDES coefficient
            30                        // Reference velocity scale
            2                         // Reference length scale
            0                         // Minimum sigma limit (0-1)
            1                         // Maximum sigma limit (0-1)
            1.0e-03                   // Limiter of B function, typically 1e-03
            1.0;                      // nut limiter (if > 1, GAM extension is active)
        .
        .
    }
    \endverbatim

Notes
  - Scheme 1 should be linear (or other low-dissipative schemes) which will
    be used in the detached/vortex shedding regions.
  - Scheme 2 should be an upwind/deferred correction/TVD scheme which will
    be used in the free-stream/Euler/boundary layer regions.
  - The scheme is compiled into a separate library, and not available to
    solvers by default.  In order to use the scheme, add the library as a
    run-time loaded library in the \$FOAM\_CASE/system/controlDict
    dictionary, e.g.:
    \verbatim
        libs (turbulenceModelSchemes);
    \endverbatim

SourceFiles
    DEShybrid.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_DEShybrid_H
#define Foam_DEShybrid_H

#include "surfaceInterpolationScheme.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "blendedSchemeBase.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class DEShybrid Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DEShybrid
:
    public surfaceInterpolationScheme<Type>,
    public blendedSchemeBase<Type>
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    // Private Data

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type>> tScheme2_;

        //- Name of the LES delta field
        word deltaName_;

        //- DES coefficient
        scalar CDES_;

        //- Reference velocity scale [m/s]
        dimensionedScalar U0_;

        //- Reference length scale [m]
        dimensionedScalar L0_;

        //- Minimum bound for sigma (0 <= sigmaMin <= 1)
        scalar sigmaMin_;

        //- Maximum bound for sigma (0 <= sigmaMax <= 1)
        scalar sigmaMax_;

        //- Limiter of B function
        scalar OmegaLim_;

        //- Limiter for modified GAM behaviour
        scalar nutLim_;

        //- Scheme constants
        scalar CH1_;
        scalar CH2_;
        scalar CH3_;
        scalar Cs_;

        //- No copy construct
        DEShybrid(const DEShybrid&) = delete;

        //- No copy assignment
        void operator=(const DEShybrid&) = delete;


    // Private Member Functions

        //- Check the scheme coefficients
        void checkValues()
        {
            if (U0_.value() <= 0)
            {
                FatalErrorInFunction
                    << "U0 coefficient must be > 0. "
                    << "Current value: " << U0_ << exit(FatalError);
            }
            if (L0_.value() <= 0)
            {
                FatalErrorInFunction
                    << "L0 coefficient must be > 0. "
                    << "Current value: " << L0_ << exit(FatalError);
            }
            if (sigmaMin_ < 0)
            {
                FatalErrorInFunction
                    << "sigmaMin coefficient must be >= 0. "
                    << "Current value: " << sigmaMin_ << exit(FatalError);
            }
            if (sigmaMax_ < 0)
            {
                FatalErrorInFunction
                    << "sigmaMax coefficient must be >= 0. "
                    << "Current value: " << sigmaMax_ << exit(FatalError);
            }
            if (sigmaMin_ > 1)
            {
                FatalErrorInFunction
                    << "sigmaMin coefficient must be <= 1. "
                    << "Current value: " << sigmaMin_ << exit(FatalError);
            }
            if (sigmaMax_ > 1)
            {
                FatalErrorInFunction
                    << "sigmaMax coefficient must be <= 1. "
                    << "Current value: " << sigmaMax_ << exit(FatalError);
            }

            if (debug)
            {
                Info<< type() << "coefficients:" << nl
                    << "    delta : " << deltaName_ << nl
                    << "    CDES : " << CDES_ << nl
                    << "    U0 : " << U0_.value() << nl
                    << "    L0 : " << L0_.value() << nl
                    << "    sigmaMin : " << sigmaMin_ << nl
                    << "    sigmaMax : " << sigmaMax_ << nl
                    << "    OmegaLim : " << OmegaLim_ << nl
                    << "    nutLim : " << nutLim_ << nl
                    << "    CH1 : " << CH1_ << nl
                    << "    CH2 : " << CH2_ << nl
                    << "    CH3 : " << CH3_ << nl
                    << "    Cs : " << Cs_ << nl
                    << endl;
            }
        }


        //- Calculate the blending factor
        tmp<surfaceScalarField> calcBlendingFactor
        (
            const VolFieldType& vf,
            const volScalarField& nut,
            const volScalarField& nu,
            const volVectorField& U,
            const volScalarField& delta
        ) const
        {
            tmp<volTensorField> tgradU = fvc::grad(U);
            const volTensorField& gradU = tgradU.cref();
            const volScalarField S(sqrt(2.0)*mag(symm(gradU)));
            const volScalarField Omega(sqrt(2.0)*mag(skew(tgradU)));
            const dimensionedScalar tau0_ = L0_/U0_;

            tmp<volScalarField> tB =
                CH3_*Omega*max(S, Omega)
               /max(0.5*(sqr(S) + sqr(Omega)), sqr(OmegaLim_/tau0_));

            tmp<volScalarField> tg = tanh(pow4(tB));

            tmp<volScalarField> tK =
                max(Foam::sqrt(0.5*(sqr(S) + sqr(Omega))), 0.1/tau0_);

            tmp<volScalarField> tlTurb =
                Foam::sqrt
                (
                    max
                    (
                        (max(nut, min(sqr(Cs_*delta)*S, nutLim_*nut)) + nu)
                       /(pow(0.09, 1.5)*tK),
                        dimensionedScalar(sqr(dimLength), Zero)
                    )
                );

            const volScalarField A
            (
                CH2_*max
                (
                    scalar(0),
                    CDES_*delta/max(tlTurb*tg, SMALL*L0_) - 0.5
                )
            );


            const word factorName(IOobject::scopedName(typeName, "Factor"));
            const fvMesh& mesh = this->mesh();

            IOobject factorIO
            (
                factorName,
                mesh.time().timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            );

            if (blendedSchemeBaseName::debug)
            {
                auto* factorPtr = mesh.getObjectPtr<volScalarField>(factorName);

                if (!factorPtr)
                {
                    factorPtr =
                        new volScalarField
                        (
                            factorIO,
                            mesh,
                            dimensionedScalar(dimless, Zero)
                        );

                    regIOobject::store(factorPtr);
                }

                auto& factor = *factorPtr;

                factor = max(sigmaMax_*tanh(pow(A, CH1_)), sigmaMin_);

                return tmp<surfaceScalarField>::New
                (
                    vf.name() + "BlendingFactor",
                    fvc::interpolate(factor)
                );
            }
            else
            {
                factorIO.registerObject(IOobjectOption::NO_REGISTER);

                volScalarField factor
                (
                    factorIO,
                    max(sigmaMax_*tanh(pow(A, CH1_)), sigmaMin_)
                );

                return tmp<surfaceScalarField>::New
                (
                    vf.name() + "BlendingFactor",
                    fvc::interpolate(factor)
                );
            }
        }


public:

    //- Runtime type information
    TypeName("DEShybrid");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        DEShybrid(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_(surfaceInterpolationScheme<Type>::New(mesh, is)),
            tScheme2_(surfaceInterpolationScheme<Type>::New(mesh, is)),
            deltaName_(is),
            CDES_(readScalar(is)),
            U0_("U0", dimLength/dimTime, readScalar(is)),
            L0_("L0", dimLength, readScalar(is)),
            sigmaMin_(readScalar(is)),
            sigmaMax_(readScalar(is)),
            OmegaLim_(readScalar(is)),
            nutLim_(readScalarOrDefault(is, scalar(1))),
            CH1_(3.0),
            CH2_(1.0),
            CH3_(2.0),
            Cs_(0.18)
        {
            checkValues();
        }

        //- Construct from mesh, faceFlux and Istream
        DEShybrid
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            deltaName_(is),
            CDES_(readScalar(is)),
            U0_("U0", dimLength/dimTime, readScalar(is)),
            L0_("L0", dimLength, readScalar(is)),
            sigmaMin_(readScalar(is)),
            sigmaMax_(readScalar(is)),
            OmegaLim_(readScalar(is)),
            nutLim_(readScalarOrDefault(is, scalar(1))),
            CH1_(3.0),
            CH2_(1.0),
            CH3_(2.0),
            Cs_(0.18)
        {
            checkValues();
        }


    // Member Functions

        //- Return the face-based blending factor
        virtual tmp<surfaceScalarField> blendingFactor
        (
             const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

            // Retrieve LES delta from the mesh database
            const auto& delta =
                mesh.lookupObject<const volScalarField>(deltaName_);

            // Retrieve turbulence model from the mesh database
            const auto* modelPtr =
                mesh.cfindObject<turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            if (modelPtr)
            {
                const auto& model = *modelPtr;

                return calcBlendingFactor
                (
                    vf, model.nut(), model.nu(), model.U(), delta
                );
            }

            FatalErrorInFunction
                << "Scheme requires a turbulence model to be present. "
                << "Unable to retrieve turbulence model from the mesh "
                << "database" << exit(FatalError);

            return nullptr;
        }


        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights(const VolFieldType& vf) const
        {
            const surfaceScalarField bf(blendingFactor(vf));

            return
                (scalar(1) - bf)*tScheme1_().weights(vf)
              + bf*tScheme2_().weights(vf);
        }


        //- Return the face-interpolate of the given cell field
        //- with explicit correction
        tmp<SurfaceFieldType> interpolate(const VolFieldType& vf) const
        {
            const surfaceScalarField bf(blendingFactor(vf));

            return
                (scalar(1) - bf)*tScheme1_().interpolate(vf)
              + bf*tScheme2_().interpolate(vf);
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme1_().corrected() || tScheme2_().corrected();
        }


        //- Return the explicit correction to the face-interpolate
        //- for the given field
        virtual tmp<SurfaceFieldType> correction(const VolFieldType& vf) const
        {
            const surfaceScalarField bf(blendingFactor(vf));

            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        (scalar(1) - bf)
                      * tScheme1_().correction(vf)
                      + bf
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        (scalar(1) - bf)
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return (bf*tScheme2_().correction(vf));
            }

            return nullptr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
