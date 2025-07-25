/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::liquid

Description
    Generic thermophysical properties class for a liquid in which the
    functions and coefficients for each property are run-time selected.

SourceFiles
    liquid.C
    liquidI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_liquid_H
#define Foam_liquid_H

#include "liquidProperties.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class liquid Declaration
\*---------------------------------------------------------------------------*/

class liquid
:
    public liquidProperties
{
    // Private Data

        autoPtr<Function1<scalar>> rho_;
        autoPtr<Function1<scalar>> pv_;
        autoPtr<Function1<scalar>> hl_;
        autoPtr<Function1<scalar>> Cp_;
        autoPtr<Function1<scalar>> h_;
        autoPtr<Function1<scalar>> Cpg_;
        autoPtr<Function1<scalar>> B_;
        autoPtr<Function1<scalar>> mu_;
        autoPtr<Function1<scalar>> mug_;
        autoPtr<Function1<scalar>> kappa_;
        autoPtr<Function1<scalar>> kappag_;
        autoPtr<Function1<scalar>> sigma_;
        autoPtr<Function1<scalar>> D_;


public:

    friend class liquidProperties;

    //- Runtime type information
    TypeName("liquid");


    // Constructors

        //- Construct from dictionary
        explicit liquid(const dictionary& dict);

        //- Copy construct
        liquid(const liquid& rhs);


        //- Return a clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return liquidProperties::Clone(*this);
        }


    // Member Functions

        //- Liquid density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Vapour pressure [Pa]
        inline scalar pv(scalar p, scalar T) const;

        //- Heat of vapourisation [J/kg]
        inline scalar hl(scalar p, scalar T) const;

        //- Liquid heat capacity [J/(kg K)]
        inline scalar Cp(scalar p, scalar T) const;

        //- Liquid Enthalpy [J/(kg)]
        inline scalar h(scalar p, scalar T) const;

        //- Ideal gas heat capacity [J/(kg K)]
        inline scalar Cpg(scalar p, scalar T) const;

        //- Second Virial Coefficient [m^3/kg]
        inline scalar B(scalar p, scalar T) const;

        //- Liquid viscosity [Pa s]
        inline scalar mu(scalar p, scalar T) const;

        //- Vapour viscosity [Pa s]
        inline scalar mug(scalar p, scalar T) const;

        //- Liquid thermal conductivity  [W/(m K)]
        inline scalar kappa(scalar p, scalar T) const;

        //- Vapour thermal conductivity  [W/(m K)]
        inline scalar kappag(scalar p, scalar T) const;

        //- Surface tension [N/m]
        inline scalar sigma(scalar p, scalar T) const;

        //- Vapour diffusivity [m2/s]
        inline scalar D(scalar p, scalar T) const;

        //- Vapour diffusivity [m2/s] with specified binary pair
        inline scalar D(scalar p, scalar T, scalar Wb) const;


    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liquidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
