/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd
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
    Foam::regionModels::pyrolysisModels::thermo

Description
    Pyrolysis model which solves only the energy equation in the region.

SourceFiles
    thermo.C

\*---------------------------------------------------------------------------*/

#ifndef regionModels_pyrolysisModels_thermo_H
#define regionModels_pyrolysisModels_thermo_H

#include "pyrolysisModel.H"
#include "volFieldsFwd.H"
#include "basicSolidChemistryModel.H"
#include "radiationModel.H"
#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

/*---------------------------------------------------------------------------*\
                          Class thermo Declaration
\*---------------------------------------------------------------------------*/

class thermo
:
    public pyrolysisModel
{
private:

    // Private member functions

        //- No copy construct
        thermo(const thermo&) = delete;

        //- No copy assignment
        void operator=(const thermo&) = delete;


protected:

    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);

        //- Read control options
        void readControls();


    // Protected data

        //- Pointer to the solid chemistry model
        autoPtr<solidThermo> solidThermo_;

        //- Pointer to radiation model
        autoPtr<radiation::radiationModel> radiation_;


        // Solution parameters

            //- Number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Maximum diffusivity
            scalar maxDiff_;


public:

    //- Runtime type information
    TypeName("thermo");


    // Constructors

        //- Construct from type name and mesh
        thermo
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& regionType
        );

        //- Construct from type name and mesh and dict
        thermo
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& regionType
        );


    //- Destructor
    virtual ~thermo();


    // Member Functions

        // Fields

            //- Return density [kg/m3]
            virtual const volScalarField& rho() const;

            //- Return const temperature [K]
            virtual const volScalarField& T() const;

            //- Return specific heat capacity [J/kg/K]
            virtual const tmp<volScalarField> Cp() const;

            //- Return the region absorptivity [1/m]
            virtual tmp<volScalarField> kappaRad() const;

            //- Return the region thermal conductivity [W/m/k]
            virtual tmp<volScalarField> kappa() const;

            //- Return the total gas mass flux to primary region [kg/m2/s]
            virtual const surfaceScalarField& phiGas() const;


        // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();

            //- Evolve the pyrolysis equations
            virtual void evolveRegion();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
