/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "Random.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Random::Random(const label seedValue)
:
    seed_(seedValue),
    generator_(seed_),
    uniform01_(),
    hasGaussSample_(false),
    gaussSample_(0)
{}


Foam::Random::Random(const Random& rnd, const bool reset)
:
    Random(rnd)
{
    if (reset)
    {
        hasGaussSample_ = false;
        gaussSample_ = 0;
        generator_.seed(seed_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::scalar Foam::Random::sample01()
{
    return scalar01();
}


template<>
Foam::label Foam::Random::sample01()
{
    return round(scalar01());
}


template<>
Foam::scalar Foam::Random::GaussNormal()
{
    if (hasGaussSample_)
    {
        hasGaussSample_ = false;
        return gaussSample_;
    }

    // Gaussian random number as per Knuth/Marsaglia.
    // Input: two uniform random numbers, output: two Gaussian random numbers.
    // cache one of the values for the next call.
    scalar rsq, v1, v2;
    do
    {
        v1 = 2*scalar01() - 1;
        v2 = 2*scalar01() - 1;
        rsq = sqr(v1) + sqr(v2);
    } while (rsq >= 1 || rsq == 0);

    const scalar fac = sqrt(-2*log(rsq)/rsq);

    gaussSample_ = v1*fac;
    hasGaussSample_ = true;

    return v2*fac;
}


template<>
Foam::label Foam::Random::GaussNormal()
{
    return round(GaussNormal<scalar>());
}


template<>
Foam::scalar Foam::Random::position
(
    const scalar& start,
    const scalar& end
)
{
    return start + scalar01()*(end - start);
}


template<>
Foam::label Foam::Random::position(const label& start, const label& end)
{
    #ifdef FULLDEBUG
    if (start > end)
    {
        FatalErrorInFunction
            << "start index " << start << " > end index " << end << nl
            << abort(FatalError);
    }
    #endif

    // Extend the upper sampling range by 1 and floor the result.
    // Since the range is non-negative, can use integer truncation
    // instead using floor().

    const label val = start + label(scalar01()*(end - start + 1));

    // Rare case when scalar01() returns exactly 1.000 and the truncated
    // value would be out of range.
    return min(val, end);
}


template<>
Foam::scalar Foam::Random::globalSample01()
{
    scalar value(0);

    if (Pstream::master())
    {
        value = scalar01();
    }

    Pstream::broadcast(value);

    return value;
}


template<>
Foam::label Foam::Random::globalSample01()
{
    label value(0);

    if (Pstream::master())
    {
        value = round(scalar01());
    }

    Pstream::broadcast(value);

    return value;
}


template<>
Foam::scalar Foam::Random::globalGaussNormal()
{
    scalar value(0);

    if (Pstream::master())
    {
        value = GaussNormal<scalar>();
    }

    Pstream::broadcast(value);

    return value;
}


template<>
Foam::label Foam::Random::globalGaussNormal()
{
    label value(0);

    if (Pstream::master())
    {
        value = GaussNormal<label>();
    }

    Pstream::broadcast(value);

    return value;
}


template<>
Foam::scalar Foam::Random::globalPosition
(
    const scalar& start,
    const scalar& end
)
{
    scalar value(0);

    if (Pstream::master())
    {
        value = position<scalar>(start, end);
    }

    Pstream::broadcast(value);

    return value;
}


template<>
Foam::label Foam::Random::globalPosition
(
    const label& start,
    const label& end
)
{
    label value(0);

    if (Pstream::master())
    {
        value = position<label>(start, end);
    }

    Pstream::broadcast(value);

    return value;
}


// ************************************************************************* //
