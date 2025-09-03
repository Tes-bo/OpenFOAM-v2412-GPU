/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

InNamespace
    Foam::IOmanip

Description
    Istream and Ostream manipulators taking arguments.

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOmanip_H
#define Foam_IOmanip_H

#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class T> class Smanip;
template<class T> class Imanip;
template<class T> class Omanip;

template<class T>
inline Istream& operator>>(Istream& is, const Smanip<T>& m);

template<class T>
inline Ostream& operator<<(Ostream& os, const Smanip<T>& m);

template<class T>
inline Istream& operator>>(Istream& is, const Imanip<T>& m);

template<class T>
inline Ostream& operator<<(Ostream& os, const Omanip<T>& m);


/*---------------------------------------------------------------------------*\
                        Class Smanip Declaration
\*---------------------------------------------------------------------------*/

//- An IOstream manipulator taking arguments
template<class T>
class Smanip
{
    T (IOstream::*_fPtr)(const T);
    T _i;

public:

    Smanip(T (IOstream::*fPtr)(const T), const T i)
    :
        _fPtr(fPtr),
        _i(i)
    {}

    friend Istream& operator>> <T>(Istream& is, const Smanip<T>& m);
    friend Ostream& operator<< <T>(Ostream& os, const Smanip<T>& m);
};


template<class T>
inline Istream& operator>>(Istream& is, const Smanip<T>& m)
{
    (is.*m._fPtr)(m._i);
    return is;
}


template<class T>
inline Ostream& operator<<(Ostream& os, const Smanip<T>& m)
{
    (os.*m._fPtr)(m._i);
    return os;
}


/*---------------------------------------------------------------------------*\
                        Class Imanip Declaration
\*---------------------------------------------------------------------------*/

//- An Istream manipulator taking arguments
template<class T>
class Imanip
{
    T (Istream::*_fPtr)(const T);
    T _i;

public:

    Imanip(T (Istream::*fPtr)(const T), const T i)
    :
        _fPtr(fPtr),
        _i(i)
    {}

    friend Istream& operator>> <T>(Istream& is, const Imanip<T>& m);
};


template<class T>
inline Istream& operator>>(Istream& is, const Imanip<T>& m)
{
    (is.*m._fPtr)(m._i);
    return is;
}


/*---------------------------------------------------------------------------*\
                        Class Omanip Declaration
\*---------------------------------------------------------------------------*/

//- An Ostream manipulator taking arguments
template<class T>
class Omanip
{
    T (Ostream::*_fPtr)(const T);
    T _i;

public:

    Omanip(T (Ostream::*fPtr)(const T), const T i)
    :
        _fPtr(fPtr),
        _i(i)
    {}

    friend Ostream& operator<< <T>(Ostream& os, const Omanip<T>& m);
};


template<class T>
inline Ostream& operator<<(Ostream& os, const Omanip<T>& m)
{
    (os.*m._fPtr)(m._i);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Smanip<std::ios_base::fmtflags> setf(std::ios_base::fmtflags flags)
{
    return Smanip<std::ios_base::fmtflags>(&IOstream::setf, flags);
}


inline Omanip<char> setfill(char fillch)
{
    return Omanip<char>(&Ostream::fill, fillch);
}


inline Omanip<IOstreamOption::streamFormat> setformat
(
    const IOstreamOption::streamFormat fmt
)
{
    return Omanip<IOstreamOption::streamFormat>(&IOstreamOption::format, fmt);
}


inline Omanip<IOstreamOption::versionNumber> setversion
(
    const IOstreamOption::versionNumber ver
)
{
    return Omanip<IOstreamOption::versionNumber>(&IOstreamOption::version, ver);
}


inline Omanip<int> setw(const int i)
{
    return Omanip<int>(&Ostream::width, i);
}


inline Omanip<int> setprecision(const int i)
{
    return Omanip<int>(&Ostream::precision, i);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
