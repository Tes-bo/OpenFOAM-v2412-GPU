/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::regionFaModels::thermalShell

Description
    Thermal-shell finite-area model. It solves the energy
    equation in 2D. The coupling with the 3D region is done through
    the \c temperatureCoupledBase, plus \c faOption is available to
    add extra sources on the shell such as \c externalHeatFluxSource etc.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        thermalShellModel   thermalShell;
        thermo
        {
            // subdictionary entries
        }

        // Optional entries
        qr                  <word>;
        thickness           <scalar>;

        // Inherited entries
        ...
        nNonOrthCorr        <int>;    // read from another dict
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                      | Type  | Reqd   | Deflt
      thermalShellModel | Type name: thermalShell   | word  | yes    | -
      thermo     | Solid thermal properties         | dictionary | yes | -
      qr         | Name of radiative heat flux field | word | no     | none
      thickness  | Uniform film thickness [m]       | scalar | choice | -
    \endtable

    The inherited entries are elaborated in:
      - \link thermalShellModel.H \endlink

SourceFiles
    thermalShell.C
    thermalShellI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionModels_thermalShell_H
#define Foam_regionModels_thermalShell_H

#include "volFieldsFwd.H"
#include "thermalShellModel.H"
#include "solidProperties.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalShell Declaration
\*---------------------------------------------------------------------------*/

class thermalShell
:
    public thermalShellModel
{
    // Private Member Functions

        //- Initialize thermalShell
        bool init(const dictionary& dict);

        //- Return radiative heat flux
        tmp<areaScalarField> qr();


protected:

    // Protected Data

        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Thermo properties

            //- Solid properties
            solidProperties thermo_;


        // Source term fields

            //- External surface energy source [J/m2/s]
            areaScalarField qs_;

            //- Film thickness [m]
            areaScalarField h_;

            //- Name of the primary region radiative flux
            const word qrName_;

            //- Uniform film thickness [m]
            scalar thickness_;


    // Protected Member Functions

        // Equations

            //- Solve energy equation
            void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermalShell");


    // Constructors

        //- Construct from components and dict
        thermalShell
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        thermalShell(const thermalShell&) = delete;

        //- No copy assignment
        void operator=(const thermalShell&) = delete;


    //- Destructor
    virtual ~thermalShell() = default;


    // Member Functions

        // Fields

            //- Return the film specific heat capacity [J/kg/K]
            const tmp<areaScalarField> Cp() const;

            //- Return density [kg/m3]
            const tmp<areaScalarField> rho() const;

            //- Return thermal conductivity [W/m/K]
            const tmp<areaScalarField> kappa() const;


        // Evolution

            //- Pre-evolve  thermal baffle
            virtual void preEvolveRegion();

            //- Evolve the thermal baffle
            virtual void evolveRegion();


       // IO

            //- Provide some feedback
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
