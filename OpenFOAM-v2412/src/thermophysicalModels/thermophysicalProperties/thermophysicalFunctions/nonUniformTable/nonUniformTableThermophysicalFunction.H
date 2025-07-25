/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::thermophysicalFunctions::nonUniformTable

Description
    Non-uniform tabulated property function that linearly interpolates between
    the values.

    To speed-up the search of the non-uniform table a uniform jump-table is
    created on construction which is used for fast indirect addressing into
    the table.

Usage
    \table
        Property    | Description
        values      | List of (temperature property) value pairs
    \endTable

    Example for the density of water between 280 and 350K
    \verbatim
    rho
    {
        type    nonUniformTable;

        values
        (
            (280 999.87)
            (300 995.1)
            (350 973.7)
        );
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef nonUniformTableThermophysicalFunction_H
#define nonUniformTableThermophysicalFunction_H

#include "thermophysicalFunction.H"
#include "MinMax.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nonUniformTable Declaration
\*---------------------------------------------------------------------------*/

class nonUniformTable
:
    public thermophysicalFunction
{
    // Private Data

        //- Name of dictionary from which this function is instantiated
        word name_;

        //- Table values
        List<Tuple2<scalar, scalar>> values_;

        //- Lowest/highest temperature in the table
        MinMax<scalar> Trange_;

        //- Temperature increment, based on Trange_ and values_.size()
        scalar deltaT_;

        //- Lookup indices into values_ table based on deltaT_
        List<label> jumpTable_;


protected:

        //- Return the lower index of the interval in the table
        //- corresponding to the given temperature.
        //  Checks the temperature range
        inline label index(scalar p, scalar T) const;


public:

    //- Runtime type information
    TypeName("nonUniformTable");


    // Constructors

        //- Construct from entry name and dictionary
        nonUniformTable(const word& name, const dictionary& dict);

        //- Construct from dictionary, using "values" for the lookup
        explicit nonUniformTable(const dictionary& dict);


    // Member Functions

        //- Return the non-uniform table of values
        const List<Tuple2<scalar, scalar>>& values() const
        {
            return values_;
        }

        //- Evaluate the function and return the result
        scalar f(scalar p, scalar T) const;

        //- Evaluate the derivative of the function and return the result
        scalar dfdT(scalar p, scalar T) const;

        //- Write the function coefficients
        void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nonUniformTableThermophysicalFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
