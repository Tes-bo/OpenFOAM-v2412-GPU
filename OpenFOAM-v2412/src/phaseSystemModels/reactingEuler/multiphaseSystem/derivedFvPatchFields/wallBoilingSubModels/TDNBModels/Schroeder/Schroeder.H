/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::TDNBModels::Schroeder

Description
    A model for departure from nucleation boiling based on
    Schroeder-Richter and Bartsch (1994) for boiling flows.

    Reference:
    \verbatim
        Schroeder-Richter, D., & Bartsch, G. (1994).
        Analytical calculation of DNB-superheating by a postulated
        thermo-mechanical effect of nucleate boiling.
        International journal of multiphase flow, 20(6), 1143-1167.
        DOI:10.1016/0301-9322(94)90060-4

        Theler, G., & Freis, D. (2011).
        Theoretical critical heat flux prediction based on non-equilibrium
        thermodynamics considerations of the subcooled boiling phenomenon.
        Mecánica Computacional, 30(19), 1713-1732.
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    TDNBModel
    {
        // Mandatory entries
        type            Schroeder;

        // Optional entries
        kg              <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                       | Type | Reqd | Deflt
      type       | Type name: Schroeder              | word | yes  | -
      kg  | Isoentropic expansion factor for ideal gases | scalar | no  | 1.666
    \endtable

SourceFiles
    Schroeder.C

\*---------------------------------------------------------------------------*/

#ifndef Schroeder_H
#define Schroeder_H

#include "TDNBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{

/*---------------------------------------------------------------------------*\
                         Class Schroeder Declaration
\*---------------------------------------------------------------------------*/

class Schroeder
:
    public TDNBModel
{
    // Private Data

        //- Isoentropic expansion factor for ideal gases
        // 5/3 monoatomic
        // 7/5 diatomic
        scalar kg_;


    // Private Member Functions

        //- No copy construct
        Schroeder(const Schroeder&) = delete;

        //- No copy assignment
        void operator=(const Schroeder&) = delete;


public:

    //- Runtime type information
    TypeName("Schroeder");


    // Constructors

        //- Construct from a dictionary
        Schroeder(const dictionary& dict);


    //- Destructor
    virtual ~Schroeder() = default;


    // Member Functions

        //- Calculate and return the departure from nulceate boiling correlation
        virtual tmp<scalarField> TDNB
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TDNBModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
