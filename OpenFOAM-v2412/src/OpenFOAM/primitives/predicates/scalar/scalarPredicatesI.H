/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::predicates::scalars::unary Foam::predicates::scalars::operation
(
    const word& opName,
    const scalar opVal,
    const scalar tol
)
{
    return operation(opNames[opName], opVal, tol);
}


inline Foam::predicates::scalars::unary Foam::predicates::scalars::operation
(
    const Tuple2<word, scalar>& op,
    const scalar tol
)
{
    return operation(opNames[op.first()], op.second(), tol);
}


inline Foam::predicates::scalars::unary Foam::predicates::scalars::operation
(
    const std::pair<word, scalar>& op,
    const scalar tol
)
{
    return operation(opNames[op.first], op.second, tol);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::predicates::scalars::contains
(
    const scalar value,
    label pos
) const
{
    return (this->find(value, pos) >= 0);
}


inline bool Foam::predicates::scalars::match(const scalar value) const
{
    return matchAny(value);
}


inline bool Foam::predicates::scalars::matchAny(const scalar value) const
{
    for (const unary& pred : *this)
    {
        if (pred(value))
        {
            return true;
        }
    }

    return false;
}


inline bool Foam::predicates::scalars::matchAll(const scalar value) const
{
    for (const unary& pred : *this)
    {
        if (!pred(value))
        {
            return false;
        }
    }

    return (!this->empty());
}


inline Foam::labelList Foam::predicates::scalars::matching
(
    const scalar value
) const
{
    labelList indices(this->size());

    label i = 0, count = 0;
    for (const unary& pred : *this)
    {
        if (pred(value))
        {
            indices[count] = i;
            ++count;
        }
        ++i;
    }
    indices.resize(count);

    return indices;
}


inline Foam::labelList Foam::predicates::scalars::matching
(
    const UList<scalar>& input,
    const bool invert
) const
{
    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (matchAny(input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


// ************************************************************************* //
