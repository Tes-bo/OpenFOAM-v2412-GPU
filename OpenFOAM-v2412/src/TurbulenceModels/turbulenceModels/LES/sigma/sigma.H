/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::LESModels::sigma

Group
    grpLESTurbulence

Description
    The sigma SGS model.

    Reference:
    \verbatim
        Nicoud, F., Toda, H. B., Cabrit, O., Bose, S., & Lee, J. (2011).
        Using singular values to build a subgrid-scale
        model for large eddy simulations.
        Physics of fluids, 23(8), 085106.
        DOI:10.1063/1.3623274
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        sigmaCoeffs
        {
            ck          0.094;
            Csigma      1.68;
        }
    \endverbatim

Note
    The default value of the Csigma constant implemented was calibrated for
    OpenFOAM using decaying isotropic turbulence and is slightly higher than
    the value suggested in the reference publication.

SourceFiles
    sigma.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LESModels_sigma_H
#define Foam_LESModels_sigma_H

#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class sigma Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class sigma
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
    // Private Member Functions

        //- No copy construct
        sigma(const sigma&) = delete;

        //- No copy assignment
        void operator=(const sigma&) = delete;


protected:

    // Protected Data

        dimensionedScalar Ck_;
        dimensionedScalar Cw_;
        dimensionedScalar Csigma_;


    // Protected Member Functions

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("sigma");


    // Constructors

        //- Construct from components
        sigma
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
    virtual ~sigma()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return SGS disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sigma.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
