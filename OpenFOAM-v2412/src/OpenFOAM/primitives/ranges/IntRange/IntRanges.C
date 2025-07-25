/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "token.H"
#include "List.H"
#include "Istream.H"
#include "Ostream.H"
#include <numeric>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
inline static Istream& input(Istream& is, IntRange<T>& range)
{
    is.readBegin("IntRange");
    is >> range.start() >> range.size();
    is.readEnd("IntRange");

    is.check(FUNCTION_NAME);
    return is;
}

template<class T>
inline static Ostream& output(Ostream& os, const IntRange<T>& range)
{
    os  << token::BEGIN_LIST
        << range.start() << token::SPACE << range.size()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::List<Foam::label> Foam::identity(const IntRange<int32_t>& range)
{
    List<label> result;

    if (range.size() > 0)
    {
        result.resize(range.size());
        std::iota(result.begin(), result.end(), range.start());
    }

    return result;
}


#if defined(WM_LABEL_SIZE) && (WM_LABEL_SIZE >= 64)
Foam::List<Foam::label> Foam::identity(const IntRange<int64_t>& range)
{
    List<label> result;

    if (range.size() > 0)
    {
        result.resize(range.size());
        std::iota(result.begin(), result.end(), range.start());
    }

    return result;
}
#endif


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, IntRange<int32_t>& range)
{
    return input(is, range);
}


Foam::Istream& Foam::operator>>(Istream& is, IntRange<int64_t>& range)
{
    return input(is, range);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const IntRange<int32_t>& range)
{
    return output(os, range);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const IntRange<int64_t>& range)
{
    return output(os, range);
}


// ************************************************************************* //
