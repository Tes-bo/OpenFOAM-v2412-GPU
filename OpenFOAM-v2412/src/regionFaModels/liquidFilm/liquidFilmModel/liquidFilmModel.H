/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::regionFaModels::liquidFilmModel

Description
    Thin film model.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        thermo                  <dict>;

        // Optional entries
        T0                      <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                      | Type  | Reqd   | Deflt
      thermo     | Liquid thermodynamics properties | dict  | yes    | -
      T0         | Reference temperature [K]       | scalar | choice | -
    \endtable

    The inherited entries are elaborated in:
      - \link liquidFilmBase.H \endlink
      - \link filmTurbulenceModel.H \endlink
      - \link injectionModelList.H \endlink
      - \link forceList.H \endlink

SourceFiles
    liquidFilmModel.C

\*---------------------------------------------------------------------------*/

#ifndef liquidFilmModel_H
#define liquidFilmModel_H

#include "volFieldsFwd.H"
#include "liquidFilmBase.H"
#include "filmTurbulenceModel.H"
#include "liquidMixtureProperties.H"
#include "injectionModelList.H"
#include "faCFD.H"
#include "forceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class liquidFilmModel Declaration
\*---------------------------------------------------------------------------*/

class liquidFilmModel
:
    public liquidFilmBase
{
protected:

    // Thermo properties

        //- Liquid thermo
        liquidMixtureProperties thermo_;

        //- Reference temperature
        scalar Tref_;


    // Fields

        //- Density [kg/m3]
        areaScalarField rho_;

        //- Dynamic viscosity [Pa.s]
        areaScalarField mu_;

        //- Film temperature [K]
        areaScalarField Tf_;

        //- Film heat capacity [J/K]
        areaScalarField Cp_;

        //- Surface tension [m/s^2]
        areaScalarField sigma_;

        //- Film rho*height [m.kg/m^3]
        areaScalarField hRho_;


    // Mass exchange sources

        //- Mass source
        areaScalarField rhoSp_;

        //- Momentum source
        areaVectorField USp_;

        //- Normal pressure by particles
        areaScalarField pnSp_;


    // Transfer fields

        //- Film mass for transfer to cloud
        volScalarField cloudMassTrans_;

        //- Parcel diameters originating from film to cloud
        volScalarField cloudDiameterTrans_;


    // General properties

        //- Turbulence model
        autoPtr<filmTurbulenceModel> turbulence_;


        // Sub-models

            //- Available mass for transfer via sub-models
            scalarField availableMass_;

            //- Cloud injection
            injectionModelList injection_;

            //- List of film forces
            forceList forces_;


public:

    //- Runtime type information
    TypeName("liquidFilmModel");


    // Constructors

        //- Construct from components and dict
        liquidFilmModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        liquidFilmModel(const liquidFilmModel&) = delete;

        //- No copy assignment
        void operator=(const liquidFilmModel&) = delete;


    //- Destructor
    virtual ~liquidFilmModel() = default;


    // Member Functions

        // Helpers

            //- Correct thermo
            void correctThermoFields();


        // Access

            //- Access const reference mu
            const areaScalarField& mu() const noexcept { return mu_; }

            //- Access const reference rho
            const areaScalarField& rho() const noexcept { return rho_; }

            //- Access const reference sigma
            const areaScalarField& sigma() const noexcept { return sigma_; }

            //- Access const reference Tf
            const areaScalarField& Tf() const noexcept { return Tf_; }

            //- Access const reference Cp
            const areaScalarField& Cp() const noexcept { return Cp_; }

            //- Access to thermo
            const liquidMixtureProperties& thermo() const noexcept
            {
                return thermo_;
            }

            //- Access to reference temperature
            scalar Tref() const noexcept { return Tref_; }


        // Transfer fields - to the primary region (lagragian injection)

            //- Return the film mass available for transfer to cloud
            virtual const volScalarField& cloudMassTrans() const;

            //- Return the parcel diameters originating from film to cloud
            virtual const volScalarField& cloudDiameterTrans() const;


        // Evolution

            //- Pre-evolve film
            virtual void preEvolveRegion();

            //- Post-evolve film
            virtual void postEvolveRegion();


        // I-O

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
