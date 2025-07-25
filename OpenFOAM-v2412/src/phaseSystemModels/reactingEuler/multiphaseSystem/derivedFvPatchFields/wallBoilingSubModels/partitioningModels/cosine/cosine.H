/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::wallBoilingModels::partitioningModels::cosine

Description
    A cosine model for wall heat flux partitioning based on
    Tentner et al. (2006) for boiling flows.

    Reference:
    \verbatim
        Tentner, A., Lo, S., Ioilev, A., Samigulin, M.,
        Ustinenko, V., Melnikov, V., Kozlov, V. (2006).
        Advances in computational fluid dynamics modeling
        of two phase flow in a boiling water reactor fuel assembly.
        In: Proc. Int. Conf. Nuclear Engineering ICONE-14.
        Miami, Florida, USA, July 17–20.
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    partitioningModel
    {
        // Mandatory entries
        type            cosine;
        alphaLiquid1    <scalar>;
        alphaLiquid0    <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                       | Type | Reqd | Deflt
      type          | Type name: cosine                 | word | yes  | -
      alphaLiquid1  | Model parameters for threshold liquid phase fraction <!--
                    -->                               | scalar | yes  | -
      alphaLiquid0  | Model parameters for threshold liquid phase fraction <!--
                    -->                               | scalar | yes  | -
    \endtable

Note
  - Proposed threshold liquid fractions:
      - alphaLiquid1 0.1
      - alphaLiquid0 0.05

SourceFiles
    cosine.C

\*---------------------------------------------------------------------------*/

#ifndef cosine_H
#define cosine_H

#include "partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                         Class cosine Declaration
\*---------------------------------------------------------------------------*/

class cosine
:
    public partitioningModel
{
    // Private Data

        //- Model parameters for threshold liquid phase fraction
        scalar alphaLiquid1_;

        //- Model parameters for threshold liquid phase fraction
        scalar alphaLiquid0_;


    // Private Member Functions

        //- No copy construct
        cosine(const cosine&) = delete;

        //- No copy assignment
        void operator=(const cosine&) = delete;


public:

    //- Runtime type information
    TypeName("cosine");


    // Constructors

        //- Construct from a dictionary
        cosine(const dictionary& dict);


    //- Destructor
    virtual ~cosine() = default;


    // Member Functions

        //- Calculate and return the wall heat-flux partitioning
        virtual tmp<scalarField> fLiquid(const scalarField& alphaLiquid) const;

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
