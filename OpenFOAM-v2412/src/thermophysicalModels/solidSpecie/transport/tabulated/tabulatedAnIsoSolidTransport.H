/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd
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
    Foam::tabulatedAnIsoSolidTransport

Description
    Transport properties package using \c Function1 type data
    for anisotropic thermal conductivity.

Usage
    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        // kappa      <Function1<vector>>;

        kappa         table
        (
            // T      kappa
            (  200    (2.56e-5 2e-5 2e-5) )
            (  350    (3.33e-5 1e-5 1e-5) )
            (  400    (4.72e-5 3e-5 3-e5) )
        );
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                        | Type   | Reqd | Deflt
      kappa     | Thermal conductivity  | Function1\<vector\> | yes  | -
    \endtable

SourceFiles
    tabulatedAnIsoSolidTransportI.H
    tabulatedAnIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tabulatedAnIsoSolidTransport_H
#define Foam_tabulatedAnIsoSolidTransport_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Thermo> class tabulatedAnIsoSolidTransport;

template<class Thermo>
Ostream& operator<<(Ostream&, const tabulatedAnIsoSolidTransport<Thermo>&);


/*---------------------------------------------------------------------------*\
                Class tabulatedAnIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedAnIsoSolidTransport
:
    public Thermo
{
    // Private Data

        //- Thermal conductivity data [W/m/K]
        autoPtr<Function1<vector>> kappa_;


    // Constructors

        //- Construct from components
        inline tabulatedAnIsoSolidTransport
        (
            const Thermo& t,
            const autoPtr<Function1<vector>>& kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline tabulatedAnIsoSolidTransport
        (
            const word&,
            const tabulatedAnIsoSolidTransport&
        );

        //- Construct from dictionary
        explicit tabulatedAnIsoSolidTransport(const dictionary& dict);

        //- Return a clone
        inline autoPtr<tabulatedAnIsoSolidTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedAnIsoSolidTransport> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- The instantiated type name
        static word typeName()
        {
            return "tabulatedAnIso<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = false;

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline vector Kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedAnIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedAnIsoSolidTransportI.H"

#ifdef NoRepository
    #include "tabulatedAnIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
