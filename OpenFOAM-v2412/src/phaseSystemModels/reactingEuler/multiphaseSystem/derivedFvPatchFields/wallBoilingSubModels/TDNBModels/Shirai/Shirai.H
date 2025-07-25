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
    Foam::wallBoilingModels::TDNBModels::Shirai

Description
    Temperature of departure from nulceate boiling correlation.

    References:
    \verbatim
        Shirai, Y., Tatsumoto, H., Shiotsu, M., Hata, K.,
        Kobayashi, H., Naruo, Y., & Inatani, Y. (2010).
        Boiling heat transfer from a horizontal flat
        plate in a pool of liquid hydrogen.
        Cryogenics, 50(6-7), 410-416.
        DOI:10.1016/j.cryogenics.2010.04.001
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    TDNBModel
    {
        // Mandatory entries
        type            Shirai;
        Tc              <scalar>;
        Pc              <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                       | Type | Reqd | Deflt
      type       | Type name: Shirai                 | word | yes  | -
      Tc         | Critical temperature            | scalar | yes  | -
      Pc         | Critical pressure               | scalar | yes  | -
    \endtable

Note
  - Correlation based on fiting data from Fig 11 from above references.
  - Suitable for liquid Helium, Nitrogen, Oxygen and Hydrogen.

SourceFiles
    Shirai.C

\*---------------------------------------------------------------------------*/

#ifndef Shirai_H
#define Shirai_H

#include "TDNBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace TDNBModels
{

/*---------------------------------------------------------------------------*\
                         Class Shirai Declaration
\*---------------------------------------------------------------------------*/

class Shirai
:
    public TDNBModel
{
    // Private Data

        //- Critical temperature
        scalar Tc_;

        //- Critical pressure
        scalar Pc_;


    // Private Member Functions

        //- No copy construct
        Shirai(const Shirai&) = delete;

        //- No copy assignment
        void operator=(const Shirai&) = delete;


public:

    //- Runtime type information
    TypeName("Shirai");


    // Constructors

        //- Construct from a dictionary
        Shirai(const dictionary& dict);


    //- Destructor
    virtual ~Shirai() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> TDNB
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

} // End namespace TDNBModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
