/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Namespace
    Foam::distributionModels

Description
    A namespace for various probability distribution model implementations.

Class
    Foam::distributionModel

Description
    A library of runtime-selectable doubly-truncated probability distribution
    models. Returns random samples based on given distribution parameters.

    Available distribution models include:
    - binned
    - exponential
    - fixedValue
    - general
    - multi-normal
    - normal
    - Rosin-Rammler
    - Mass-based Rosin-Rammler
    - uniform

SourceFiles
    distributionModel.C
    distributionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModel_H
#define distributionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class distributionModel Declaration
\*---------------------------------------------------------------------------*/

class distributionModel
{
protected:

    // Protected Data

        //- Coefficients dictionary
        const dictionary distributionModelDict_;

        //- Reference to the random number generator
        Random& rndGen_;

        //- Minimum of the distribution
        scalar minValue_;

        //- Maximum of the distribution
        scalar maxValue_;


    // Protected Member Functions

        //- Check that the distribution model is valid
        virtual void check() const;


public:

    //-Runtime type information
    TypeName("distributionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        distributionModel,
        dictionary,
        (
            const dictionary& dict,
            Random& rndGen
        ),
        (dict, rndGen)
    );


    // Constructors

        //- Construct from dictionary
        distributionModel
        (
            const word& name,
            const dictionary& dict,
            Random& rndGen
        );

        //- Copy construct
        distributionModel(const distributionModel& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const = 0;


    //- Selector
    static autoPtr<distributionModel> New
    (
        const dictionary& dict,
        Random& rndGen
    );


    //- Destructor
    virtual ~distributionModel() = default;


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const = 0;

        //- Return the minimum of the distribution
        virtual scalar minValue() const;

        //- Return the maximum of the distribution
        virtual scalar maxValue() const;

        //- Return the theoretical mean of the distribution, or
        //- the arithmetic mean of the distribution data
        virtual scalar meanValue() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
