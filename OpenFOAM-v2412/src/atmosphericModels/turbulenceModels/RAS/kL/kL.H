/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::RASModels::kL

Group
    grpRASTurbulence

Description
    A one-equation (turbulent kinetic energy \c k) turbulence closure
    model for incompressible and compressible geophysical applications.

    Turbulent kinetic energy (\c k) is computed with a transport equation
    and the turbulent length scale (\c L) is computed with an algebraic
    expression which depends on the local stratification.

    References:
    \verbatim
        Standard model (tag:A):
            Axell, L. B., & Liungman, O. (2001).
            A one-equation turbulence model for geophysical applications:
            comparison with data and the k−ε model.
            Environmental Fluid Mechanics, 1(1), 71-106.
            DOI:10.1023/A:1011560202388

        Canopy-related models (tag:W):
            Wilson, J. D., & Flesch, T. K. (1999).
            Wind and remnant tree sway in forest cutblocks.
            III. A windflow model to diagnose spatial variation.
            Agricultural and Forest Meteorology, 93(4), 259-282.
            DOI:10.1016/S0168-1923(98)00121-X
    \endverbatim

Usage
    Example by using \c constant/turbulenceProperties:
    \verbatim
    RAS
    {
        // Mandatory entries
        RASModel        kL;

        // Optional entries
        kLCoeffs
        {
            kappa       <scalar>;
            sigmak      <scalar>;
            beta        <scalar>;
            Cmu0        <scalar>;
            Lmax        <scalar>;
            CbStable    <scalar>;
            CbUnstable  <scalar>;
        }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                     | Type   | Reqd | Deflt
      RASModel   | Type name: kL                   | word   | yes  | -
      kappa      | von Karman constant             | scalar | no   | 0.41
      sigmak     | Empirical model coefficient     | scalar | no   | 1.0
      beta | Thermal expansion coefficient [1/K]   | scalar | no   | 3.3e-3
      Cmu0       | Empirical model coefficient     | scalar | no   | 0.556
      Lmax   | Maximum mixing-length scale [m]     | scalar | no   | GREAT
      CbStable   | Stable stratification constant  | scalar | no   | 0.25
      CbUnstable | Unstable stratification constant | scalar | no  | 0.35
    \endtable

    The inherited entries are elaborated in:
      - \link eddyViscosity.H \endlink

    \heading Input fields (mandatory)
    \plaintable
        k         | Turbulent kinetic energy                [m2/s2]
        T         | Potential temperature                   [K]
    \endplaintable

    \heading Input fields (optional)
    \plaintable
        canopyHeight    | Canopy height                     [m]
        Cd              | Plant canopy drag coefficient     [-]
        LAD             | Leaf area density                 [1/m]
        Rt              | Turbulent Richardson number       [-]
        L               | Characteristic length scale       [m]
    \endplaintable

Note
  - Optional input fields can/should be input
    by using \c readFields function object.

SourceFiles
    kL.C

\*---------------------------------------------------------------------------*/

#ifndef kL_H
#define kL_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kL Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kL
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- Return modified stability function of Launder,
        //- calibrated with the results of Hogstrom (A:Eq. 31)
        tmp<volScalarField> Cmu() const;

        //- Return modified stability function of Launder,
        //- calibrated with the results of Hogstrom (A:Eq. 32)
        tmp<volScalarField> CmuPrime() const;

        //- Return eddy diffusivity (A:Eq. 12)
        tmp<volScalarField> nutPrime() const;

        //- Return turbulent kinetic energy dissipation rate due to canopy
        tmp<volScalarField> epsilonCanopy() const;

        //- Return turbulent kinetic energy
        //- dissipation rate due to plains and canopy
        tmp<volScalarField> epsilon() const;

        //- Return canopy height
        tmp<volScalarField> canopyHeight() const;

        //- Return characteristic length scale
        tmp<volScalarField> L() const;

        //- Modify characteristic length scale
        //- according to the specified stratification
        void stratification(const volScalarField& fVB);

        // Generated Methods

            //- No copy construct
            kL(const kL&) = delete;

            //- No copy assignment
            void operator=(const kL&) = delete;


protected:

    // Protected Data

        // Model coefficients

            //- von Karman constant
            dimensionedScalar kappa_;

            //- Empirical model coefficient
            dimensionedScalar sigmak_;

            //- Thermal expansion coefficient [1/K]
            dimensionedScalar beta_;

            //- Empirical model coefficient
            dimensionedScalar Cmu0_;

            //- Maximum mixing-length scalar [m]
            dimensionedScalar Lmax_;

            //- Stable stratification constant
            dimensionedScalar CbStable_;

            //- Unstable stratification constant
            dimensionedScalar CbUnstable_;


        // Fields

            //- Turbulent kinetic energy    [m2/s2]
            volScalarField k_;

            //- Characteristic length scale [m]
            volScalarField L_;

            //- Turbulent Richardson number [-]
            volScalarField Rt_;

            //- Gravitational acceleration  [m2/s2]
            const uniformDimensionedVectorField& g_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        //- Correct the turbulence viscosity
        virtual void correctNut();

        //- Add explicit source for turbulent kinetic energy
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kL");


    // Constructors

        //- Construct from components
        kL
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kL() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kL.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
