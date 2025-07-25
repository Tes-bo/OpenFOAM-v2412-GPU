/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Print some numerical limits.

\*---------------------------------------------------------------------------*/

#include <cmath>
#include <limits>
#include "int.H"
#include "uint.H"
#include "scalar.H"
#include "word.H"
#include "IOstreams.H"

using namespace Foam;

std::ostream& print(const char* tag, float val)
{
    std::cout
        << tag << val
        << " 0x" << std::hex << *reinterpret_cast<const uint32_t*>(&val);
    return std::cout;
}


std::ostream& print(const char* tag, double val)
{
    std::cout
        << tag << val
        << " 0x" << std::hex << *reinterpret_cast<const uint64_t*>(&val);
    return std::cout;
}


// Have (float|double)Scalar(GREAT|SMALL|..)

#define PrintFloatLimits(FloatType, NanFunction)                              \
{                                                                             \
    print("max:", std::numeric_limits<FloatType>::max());                     \
    print(" VGREAT:", FloatType##ScalarVGREAT);                               \
    print(" ROOTVGREAT:", FloatType##ScalarROOTVGREAT) << nl;                 \
                                                                              \
    print("min:", std::numeric_limits<FloatType>::min());                     \
    print(" VSMALL:", FloatType##ScalarVSMALL);                               \
    print(" ROOTVSMALL:", FloatType##ScalarROOTVSMALL) << nl;                 \
                                                                              \
    print("epsilon:", std::numeric_limits<FloatType>::epsilon());             \
    print(" SMALL:", FloatType##ScalarSMALL);                                 \
    print(" ROOTSMALL:", FloatType##ScalarROOTSMALL) << nl;                   \
                                                                              \
    print("1/epsilon:", 1/std::numeric_limits<FloatType>::epsilon());         \
    print(" GREAT:", FloatType##ScalarGREAT);                                 \
    print(" ROOTGREAT:", FloatType##ScalarROOTGREAT) << nl;                   \
                                                                              \
    print("nan:", std::NanFunction("nan")) << nl;                             \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    //NONE cout<<"int16:" << pTraits<int16_t>::max << nl;
    cout<< "=max=" << nl;
    cout<< "uint8:" << pTraits<uint8_t>::max << nl;
    cout<< "int16:" << std::numeric_limits<int16_t>::max() << nl;
    cout<< "int32:" << pTraits<int32_t>::max << nl;
    cout<< "uint32:" << pTraits<uint32_t>::max << nl;
    cout<< "int64:" << pTraits<int64_t>::max << nl;
    cout<< "uint64:" << pTraits<uint64_t>::max << nl;

    cout<< nl;

    cout<< "int16:" << std::numeric_limits<int16_t>::max() << nl;
    cout<< "int32:" << pTraits<int32_t>::max << nl;
    cout<< "uint32:" << pTraits<uint32_t>::max << nl;
    cout<< "int64:" << pTraits<int64_t>::max << nl;
    cout<< "uint64:" << pTraits<uint64_t>::max << nl;

    cout<< nl << "=digits=" << nl;

    cout<< "int16:" << std::numeric_limits<int16_t>::digits << nl;
    cout<< "int32:" << std::numeric_limits<int32_t>::digits << nl;
    cout<< "uint32:" << std::numeric_limits<uint32_t>::digits << nl;
    cout<< "int64:" << std::numeric_limits<int64_t>::digits << nl;
    cout<< "uint64:" << std::numeric_limits<uint64_t>::digits << nl;
    cout<< "float:"  << std::numeric_limits<float>::digits << nl;
    cout<< "double:" << std::numeric_limits<double>::digits << nl;
    cout<< "long double:" << std::numeric_limits<long double>::digits << nl;

    // std::nanl (long double)
    cout<< nl;
    cout<< nl << "=float=" << nl;
    PrintFloatLimits(float, nanf);

    cout<< nl << "=double=" << nl;
    PrintFloatLimits(double, nan);

    cout<< "---\nEnd\n" << std::endl;

    return 0;
}


// ************************************************************************* //
