/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::distributionModels::fixedValue

Description
    Particle-size distribution model wherein samples are given fixed values.

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        injectionModels
        {
            <name>
            {
                ...

                sizeDistribution
                {
                    type        fixedValue;
                    fixedValueDistribution
                    {
                        value   <value>;
                    }
                }
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                       | Type | Reqd | Deflt
      type         | Type name: fixedValue             | word | yes  | -
      fixedValueDistribution | Distribution settings   | dict | yes  | -
      value        | Fixed value for size              | scalar | yes | -
    \endtable

SourceFiles
    fixedValue.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_fixedValue_H
#define distributionModels_fixedValue_H

#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{
/*---------------------------------------------------------------------------*\
                         Class fixedValue Declaration
\*---------------------------------------------------------------------------*/

class fixedValue
:
    public distributionModel
{
    // Private Data

        //- Fixed value for size
        scalar value_;


public:

    //- Runtime type information
    TypeName("fixedValue");


    // Constructors

        //- Construct from components
        fixedValue(const dictionary& dict, Random& rndGen);

        //- Copy construct
        fixedValue(const fixedValue& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new fixedValue(*this));
        }

        //- No copy assignment
        void operator=(const fixedValue&) = delete;


    //- Destructor
    virtual ~fixedValue() = default;


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Return the minimum of the distribution
        virtual scalar minValue() const;

        //- Return the maximum of the distribution
        virtual scalar maxValue() const;

        //- Return the theoretical mean of the distribution
        virtual scalar meanValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
