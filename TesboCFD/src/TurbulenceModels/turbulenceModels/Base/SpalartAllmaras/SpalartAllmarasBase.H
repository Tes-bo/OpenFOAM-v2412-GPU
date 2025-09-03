/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::LESModels::SpalartAllmarasBase

Group
    grpDESTurbulence

Description
    Base class to handle various characteristics for \c SpalartAllmaras based
    LES/DES turbulence models for incompressible and compressible flows.

    References:
    \verbatim
        Standard model:
            Spalart, P.R., & Allmaras, S.R. (1994).
            A one-equation turbulence model for aerodynamic flows.
            La Recherche Aerospatiale, 1, 5-21.

        Standard model:
            Spalart, P. R., Jou, W. H., Strelets, M., & Allmaras, S. R. (1997).
            Comments on the feasibility of LES for wings, and on a hybrid
            RANS/LES approach.
            Advances in DNS/LES, 1, 4-8.

        Estimation expression for k and epsilon (tag:B), Eq. 4.50:
            Bourgoin, A. (2019).
            Bathymetry induced turbulence modelling the
            Alderney Race site: regional approach with TELEMAC-LES.
            Normandie Université.

        Estimation expressions for omega (tag:P):
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531
    \endverbatim

SourceFiles
    SpalartAllmarasBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpalartAllmarasBase_H
#define Foam_SpalartAllmarasBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SpalartAllmarasBase Declaration
\*---------------------------------------------------------------------------*/

template<class BasicEddyViscosityModel>
class SpalartAllmarasBase
:
    public BasicEddyViscosityModel
{
    // Private Member Functions

        //- No copy construct
        SpalartAllmarasBase(const SpalartAllmarasBase&) = delete;

        //- No copy assignment
        void operator=(const SpalartAllmarasBase&) = delete;


protected:

    // Protected Data

        // Model constants

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cs_;
            dimensionedScalar ck_;

            Switch ft2_;
            dimensionedScalar Ct3_;
            dimensionedScalar Ct4_;


        // Fields

            //- Modified kinematic viscosity [m^2/s]
            volScalarField nuTilda_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        tmp<volScalarField> ft2(const volScalarField& chi) const;

        tmp<volScalarField> Omega(const volTensorField& gradU) const;

        tmp<volScalarField> r
        (
            const volScalarField& nur,
            const volScalarField& Stilda,
            const volScalarField& dTilda
        ) const;

        tmp<volScalarField::Internal> fw
        (
            const volScalarField& Stilda,
            const volScalarField& dTilda
        ) const;

        virtual tmp<volScalarField> Stilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU,
            const volScalarField& dTilda
        ) const;

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& chi,
            const volScalarField& fv1,
            const volTensorField& gradU
        ) const = 0;

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();


public:

    typedef typename BasicEddyViscosityModel::alphaField alphaField;
    typedef typename BasicEddyViscosityModel::rhoField rhoField;
    typedef typename BasicEddyViscosityModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        SpalartAllmarasBase
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~SpalartAllmarasBase() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return the (estimated) turbulent kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the (estimated) turbulent kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the (estimated) specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Return the modified kinematic viscosity
        tmp<volScalarField> nuTilda() const
        {
            return nuTilda_;
        }

        //- Correct nuTilda and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
