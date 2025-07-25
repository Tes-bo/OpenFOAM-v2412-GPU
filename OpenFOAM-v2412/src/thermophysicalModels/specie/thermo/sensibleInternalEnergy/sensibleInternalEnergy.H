/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::sensibleInternalEnergy

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the sensible internal energy
    functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleInternalEnergy_H
#define sensibleInternalEnergy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sensibleInternalEnergy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sensibleInternalEnergy
{

public:

    //- Constructor
    sensibleInternalEnergy()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleInternalEnergy";
        }


        // Fundamental properties

            static word energyName()
            {
                return "e";
            }

            //- Heat capacity at constant volume [J/(kg K)]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar cv = thermo.Cv(p, T);
                return cv;
                #else
                return thermo.Cv(p, T);
                #endif
            }

            //- Ratio of specific heats Cp/Cv []
            scalar CpByCpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar gamma = thermo.gamma(p, T);
                return gamma;
                #else
                return thermo.gamma(p, T);
                #endif
            }

            //- Sensible internal energy [J/kg]
            scalar HE
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar es = thermo.Es(p, T);
                return es;
                #else
                return thermo.Es(p, T);
                #endif
            }

            //- Temperature from sensible internal energy given an initial
            //- temperature T0 [K]
            scalar THE
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar tes = thermo.TEs(e, p, T0);
                return tes;
                #else
                return thermo.TEs(e, p, T0);
                #endif
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
