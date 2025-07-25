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
    Foam::solidIsothermalReactionRate

Description
    Isothermal reaction rate for solids

Usage
    Minimal example:
    \verbatim
    {
        // Mandatory entries
        C       <scalar>;
        Cp      <scalar>;
        Tpc     <scalar>;
        Elat    <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description             | Type   | Reqd | Deflt
      C         | Model constant          | scalar | yes  | -
      Cp | Specific heat of solid material (assumed constant) | scalar | yes | -
      Tpc       | Phase change temperature | scalar | yes | -
      Elat      | Latent heat of phase change | scalar | yes | -
    \endtable

Note
  - \c Elat, the latent energy is defined in the present
    implementation as a negative value (using the convention
    endothermic process has negative enthalpy.)

SourceFiles
    solidIsothermalReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_solidIsothermalReactionRate_H
#define Foam_solidIsothermalReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class solidIsothermalReactionRate;

Ostream& operator<<(Ostream&, const solidIsothermalReactionRate&);


/*---------------------------------------------------------------------------*\
                Class solidIsothermalReactionRate Declaration
\*---------------------------------------------------------------------------*/

class solidIsothermalReactionRate
{
    // Private Data

        //- Model constant
        scalar C_;

        //- Specific heat of solid material (assumed constant)
        scalar Cp_;

        //- Phase change temperature
        scalar Tpc_;

        //- Latent heat of phase change
        scalar Elat_;


public:

    // Constructors

        //- Construct from components
        inline solidIsothermalReactionRate
        (
            const scalar C,
            const scalar Cp,
            const scalar Tpc,
            const scalar Elat
        );


        //- Construct from dictionary
        inline solidIsothermalReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    virtual ~solidIsothermalReactionRate() = default;


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Isothermal";
        }

        //- Return reaction rate constant
        //  frequency of collisions resulting in a reaction
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const solidIsothermalReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidIsothermalReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
