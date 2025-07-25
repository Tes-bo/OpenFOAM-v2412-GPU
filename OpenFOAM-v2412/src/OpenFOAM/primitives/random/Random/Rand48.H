/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::Rand48

Description
    A pseudo random number generator using the linear congruential algorithm
    with the following parameters:
    \verbatim
        Xn+1 = (aXn + c) mod m, where n >= 0

        a = 0x5DEECE66D,
        c = 0xB
        m = 2**48
    \endverbatim

    It delivers results identical to the \c lrand48() function that is
    available on some systems.

\*---------------------------------------------------------------------------*/

#ifndef Rand48_H
#define Rand48_H

#include <cstdint>
#include <random>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Rand48 Declaration
\*---------------------------------------------------------------------------*/

class Rand48
{
    //- It is a linear_congruential_engine
    typedef std::linear_congruential_engine
    <
        uint64_t,
        0x5DEECE66D,
        0xB,
        (uint64_t(1) << 48)
    > engine;


    // Private Data

        //- The calculation engine
        engine engine_;


    // Private Member Functions

        //- Convert integers
        static constexpr uint64_t convert(uint32_t x) noexcept
        {
            return (static_cast<uint64_t>(x) << 16) | 0x330e;
        }

public:

    //- The type of the generated random value.
    typedef uint32_t result_type;

    //- The default seed
    static constexpr result_type default_seed = 1u;


    // Constructors

        //- Construct with specified or default seed
        explicit Rand48(uint32_t val = default_seed)
        :
            engine_(convert(val))
        {}


    // Member Functions

        //- The smallest value that the generator can produce
        static constexpr uint32_t min() { return 0; }

        //- The largest value that the generator can produce
        static constexpr uint32_t max() { return 0x7FFFFFFF; }

        //- Reseed the generator
        void seed(uint32_t val = default_seed)
        {
            engine_.seed(convert(val));
        }

        //- Advance the state of the generator by \c z notches.
        void discard(unsigned long long z)
        {
            engine_.discard(z);
        }

        //- Get the next random number in the sequence and return its raw
        //- 48-bit value
        uint64_t raw()
        {
            return engine_();
        }


    // Member Operators

        //- Get the next random number in the sequence
        uint32_t operator()()
        {
            return static_cast<uint32_t>(engine_() >> 17);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
