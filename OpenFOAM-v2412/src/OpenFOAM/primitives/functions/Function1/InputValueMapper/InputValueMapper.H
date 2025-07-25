/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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
    Foam::Function1Types::InputValueMapper

Description
    Function1 wrapper that maps the input value prior to it being used by
    another Function1.

    Example usage for limiting a polynomial:
    \verbatim
        <entryName>
        {
            type            inputValueMapper;
            mode            minMax;

            min             0.4;
            max             1.4;

            value polynomial
            (
                (5 1)
                (-2 2)
                (-2 3)
                (1 4)
            );
        }
    \endverbatim

    Here the return value will be:
    - poly(0.4) for x <= 0.4;
    - poly(1.4) for x >= 1.4; and
    - poly(x) for 0.4 < x < 1.4.


    Example usage for supplying a patch mass flux for a table lookup:
    \verbatim
        <entryName>
        {
            type            inputValueMapper;
            mode            function;

            function
            {
                type            functionObjectValue;
                functionObject  surfaceFieldValue1;
                functionObjectResult sum(outlet,phi);
            }

            value
            {
                type        table;
                file        "<system>/fanCurve.txt";
            }
        }
    \endverbatim

    Where:
    \table
        Property | Description                                  | Required
        mode     | Mapping mode (see below)                     | yes
        function | Mapping Function1                            | no*
        min      | Minimum input value                          | no*
        max      | Maximum input value                          | no*
        value    | Function of type Function1<Type>             | yes
    \endtable

    Mapping modes include
    - none     : the input value is simply passed to the 'value' Function1
    - function : the input value is passed through the 'function' Function1
                 before being passed to the 'value' Function1
    - minMax   : limits the input value to 'min' and 'max' values before being
                 passed to the 'value' Function1

Note
    Replaces the LimitRange Function1 (v2106 and earlier)

SourceFiles
    InputValueMapper.C
    InputValueMapperI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_InputValueMapper_H
#define Foam_Function1Types_InputValueMapper_H

#include "Function1.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                      Class InputValueMapper Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class InputValueMapper
:
    public Function1<Type>
{
public:

    // Public Enumerations

        //- Input value mapping mode
        enum class mappingMode
        {
            NONE,
            FUNCTION1,
            MINMAX
        };

        //- Names for the input value  mapping modes
        static const Enum<mappingMode> mappingModeNames_;


private:

    // Private Data

        //- Mapping mode
        mappingMode mappingMode_;

        // Function1 mapping

            //- Input value Function1
            autoPtr<Function1<scalar>> mappingValuePtr_;


        // Min/max mapping

            //- Minimum input value
            scalar min_;

            //- Maximum input value
            scalar max_;


        // Function being wrapped

            //- Value function
            autoPtr<Function1<Type>> value_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);


public:

    //- Runtime type information
    TypeName("inputValueMapper");


    // Generated Methods

        //- No copy assignment
        void operator=(const InputValueMapper<Type>&) = delete;


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        InputValueMapper
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit InputValueMapper(const InputValueMapper<Type>& rhs);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~InputValueMapper() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const;

        //- Integrate between two (scalar) values
        virtual inline Type integrate(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InputValueMapperI.H"

#ifdef NoRepository
    #include "InputValueMapper.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
