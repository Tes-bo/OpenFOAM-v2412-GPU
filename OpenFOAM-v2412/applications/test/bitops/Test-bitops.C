/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

Description
    Test some bit-operations.

\*---------------------------------------------------------------------------*/

#include "bool.H"
#include "BitOps.H"
#include "IOstreams.H"
#include "stdFoam.H"

#include <type_traits>
#include <limits>


namespace Foam
{
namespace Detail
{

template<typename UIntType, UIntType v, unsigned int n>
struct bitops_setlower
:
    std::integral_constant
    <
        UIntType,
        v | (v >> n) | bitops_setlower<UIntType, v | (v >> n), (n >> 1)>::value
    >
{};

template<typename UIntType, UIntType v>
struct bitops_setlower<UIntType, v, 1>
:
    std::integral_constant<UIntType, v | (v >> 1)>
{};


template<size_t N>
struct bitops_topbit
:
    std::integral_constant
    <
        size_t,
        bitops_topbit<(N >> 1)>{} + 1
    >
{};

template<> struct bitops_topbit<2> : std::integral_constant<size_t,1> {};
template<> struct bitops_topbit<1> : std::integral_constant<size_t,1> {};
template<> struct bitops_topbit<0> : std::integral_constant<size_t,0> {};


template<typename UIntType, UIntType v>
struct pow2ceil
:
    public std::integral_constant
    <
        UIntType,
        bitops_setlower
        <
            UIntType,
            v - 1,
            (std::numeric_limits<UIntType>::digits >> 1)
        >::value + 1
    >
{};

template<size_t N>
struct pow2ceil_shift
:
    std::integral_constant
    <
        size_t,
        bitops_topbit<pow2ceil<size_t, N>::value>::value
    >
{};

}
}


template<typename UIntType, UIntType v>
struct pow2ceil
:
    Foam::Detail::pow2ceil<UIntType,v>
{};


template<size_t N>
struct pow2topbit
:
    Foam::Detail::pow2ceil_shift<N>
{};


using namespace Foam;

template<size_t N>
void printTopbit()
{
    std::cout
        << "pow2ceil<" << N << "> = "
        << pow2ceil<size_t, N>::value
        << " shift = " << pow2topbit<N>::value << '\n';
}


template<class T, int Offset = 19>
void printOffset()
{
    std::cout
        << "pow2ceil of " << typeid(T).name() << " <" << sizeof(T) << "> = "
        << pow2ceil<size_t, sizeof(T)>::value
        << " shift = " << pow2topbit<sizeof(T)>::value << '\n';
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<<"pow2ceil<10>: "
        << pow2ceil<unsigned, 10>::value << nl;

    Info<<"pow2ceil<16>: "
        << pow2ceil<unsigned, 16>::value << nl;

    Info<<"pow2ceil<369>: "
        << pow2ceil<unsigned, 369>::value << nl;

    printTopbit<0>();
    printTopbit<1>();
    printTopbit<2>();
    printTopbit<3>();
    printTopbit<4>();
    printTopbit<5>();
    printTopbit<6>();
    printTopbit<7>();
    printTopbit<8>();
    printTopbit<9>();
    printTopbit<10>();
    printTopbit<11>();
    printTopbit<15>();
    printTopbit<16>();
    printTopbit<17>();
    printTopbit<18>();
    printTopbit<29>();
    printTopbit<30>();
    printTopbit<31>();
    printTopbit<32>();
    printTopbit<33>();
    printTopbit<4095>();
    printTopbit<4096>();
    printTopbit<4097>();

    printOffset<double>();


    Info<<nl << "Test repeat_value" << nl << nl;

    Info<< BitOps::bitInfo<unsigned>(BitOps::repeat_value<unsigned, 3>(1u))
        << nl;

    Info<< BitOps::bitInfo<unsigned>(BitOps::repeat_value<unsigned, 1>(1u))
        << nl;


    Info << "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
