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
    Foam::wallBoilingModels::nucleateFluxModels::exponential

Description
    Nucleate flux sub-cooling correlation

    References:
    \verbatim
        Frost W. & Dzakowic G. S. (1967)
        An extension of the methods of predicting incipient
        boiling on commercially finished surfaces.
        ASME AIChE Heat Transfer Conf 67-HT-61, Seattle.

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
    nucleateFluxModel
    {
        // Mandatory entries
        type            exponential;

        // Optional entries
        a               <scalar>;
        b               <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                       | Type | Reqd | Deflt
      type       | Type name: exponential            | word | yes  | -
      a          | Pre-factor coefficient          | scalar | no   | 6309
      b          | Exponent coefficient            | scalar | no   | 2.52
    \endtable

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "nucleateFluxModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleateFluxModels
{

/*---------------------------------------------------------------------------*\
                         Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public nucleateFluxModel
{
    // Private Data

        //- Pre-factor coefficient
        scalar a_;

        //- Exponent coefficient
        scalar b_;


    // Private Member Functions

        //- No copy construct
        exponential(const exponential&) = delete;

        //- No copy assignment
        void operator=(const exponential&) = delete;


public:

    //- Runtime type information
    TypeName("exponential");


    // Constructors

        //- Construct from a dictionary
        exponential(const dictionary& dict);


    //- Destructor
    virtual ~exponential() = default;


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
