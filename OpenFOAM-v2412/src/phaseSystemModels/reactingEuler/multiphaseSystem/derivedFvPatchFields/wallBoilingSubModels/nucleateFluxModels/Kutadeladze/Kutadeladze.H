/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::nucleateFluxModels::Kutadeladze

Description
    Nucleate flux sub-cooling correlation

    References:
    \verbatim
        Wang, L., Li, Y., Zhang, F., Xie, F., & Ma, Y. (2016).
        Correlations for calculating heat transfer of hydrogen pool boiling.
        International Journal of Hydrogen Energy, 41(38), 17118-17131.
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    nucleateFluxModel
    {
        // Mandatory entries
        type            Kutadeladze;

        // Optional entries
        Cn              <scalar>;
        an              <scalar>;
        bn              <scalar>;
        n               <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                       | Type | Reqd | Deflt
      type       | Type name: Kutadeladze            | word | yes  | -
      Cn         | Model coefficient               | scalar | no   | 5.66e-10
      an      | Coefficient for deltaT sub-cooling | scalar | no   | 2.5
      bn      | Exponent for n                     | scalar | no   | 1
      n       | Nucleating boiling paramemeter     | scalar | no   | 1
    \endtable

SourceFiles
    Kutadeladze.C

\*---------------------------------------------------------------------------*/

#ifndef Kutadeladze_H
#define Kutadeladze_H

#include "nucleateFluxModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleateFluxModels
{

/*---------------------------------------------------------------------------*\
                         Class Kutadeladze Declaration
\*---------------------------------------------------------------------------*/

class Kutadeladze
:
    public nucleateFluxModel
{
    // Private Data

        //- Model coefficient
        scalar Cn_;

        //- Coefficient for deltaT sub-cooling
        scalar an_;

        //- Exponent for n
        scalar bn_;

        //- Nucleating boiling paramemeter
        scalar n_;


    // Private Member Functions

        //- No copy construct
        Kutadeladze(const Kutadeladze&) = delete;

        //- No copy assignment
        void operator=(const Kutadeladze&) = delete;


public:

    //- Runtime type information
    TypeName("Kutadeladze");


    // Constructors

        //- Construct from a dictionary
        Kutadeladze(const dictionary& dict);


    //- Destructor
    virtual ~Kutadeladze() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> qNucleate
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleateFluxModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
