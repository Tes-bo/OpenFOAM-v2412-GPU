/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Print values of predefined dimensionSets, and some other tests

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "IOmanip.H"
#include <tuple>

using namespace Foam;

// Macros to stringify macro contents.
#define STRINGIFY(content)      #content
#define STRING_QUOTE(input)     STRINGIFY(input)

#define PRINT_DIMSET(arg)       \
    Info<< STRING_QUOTE(arg) << "  " << arg << nl


bool hadDimensionError
(
    const std::tuple<bool, dimensionSet, dimensionSet>& input,
    bool dimsOk,
    std::string errMsg
)
{
    if (dimsOk)
    {
        if (std::get<0>(input))
        {
            Info<< "(pass) dimension check ok ";
        }
        else
        {
            Info<< "(fail) unexpected success for dimension check ";
        }
        Info<< std::get<1>(input) << " == " << std::get<2>(input) << nl;
    }
    else
    {
        if (std::get<0>(input))
        {
            Info<< "(fail) unexpected";
        }
        else
        {
            Info<< "(pass) expected";
        }

        Info<< " failure" << nl << errMsg.c_str() << nl;
    }

    return (std::get<0>(input) != dimsOk);
}


unsigned checkDimensions
(
    std::initializer_list
    <
        std::tuple<bool, dimensionSet, dimensionSet>
    > tests
)
{
    Info<< nl << "Verify dimension checks" << nl << nl;

    unsigned nFail = 0;
    std::string errMsg;

    // Expect some failures
    const bool oldThrowingError = FatalError.throwing(true);

    for
    (
        const std::tuple<bool, dimensionSet, dimensionSet>& test
      : tests
    )
    {
        const bool expected = std::get<0>(test);
        const dimensionSet& a = std::get<1>(test);
        const dimensionSet& b = std::get<2>(test);

        bool dimsOk = false;

        try
        {
            min(a, b);
            dimsOk = true;
        }
        catch (const Foam::error& err)
        {
            errMsg = err.message();
        }

        if (expected != dimsOk)
        {
            ++nFail;
        }

        hadDimensionError(test, dimsOk, errMsg);
    }

    FatalError.throwing(oldThrowingError);

    return nFail;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Predefined dimensionSets" << nl << nl;

    PRINT_DIMSET(dimless);
    PRINT_DIMSET(dimMass);
    PRINT_DIMSET(dimLength);
    PRINT_DIMSET(dimTime);
    PRINT_DIMSET(dimTemperature);
    PRINT_DIMSET(dimMoles);
    PRINT_DIMSET(dimCurrent);
    PRINT_DIMSET(dimLuminousIntensity);

    PRINT_DIMSET(dimArea);
    PRINT_DIMSET(dimVolume);
    PRINT_DIMSET(dimVol);

    PRINT_DIMSET(dimVelocity);
    PRINT_DIMSET(dimAcceleration);

    PRINT_DIMSET(dimDensity);
    PRINT_DIMSET(dimForce);
    PRINT_DIMSET(dimEnergy);
    PRINT_DIMSET(dimPower);

    PRINT_DIMSET(dimPressure);
    PRINT_DIMSET(dimCompressibility);
    PRINT_DIMSET(dimGasConstant);
    PRINT_DIMSET(dimSpecificHeatCapacity);
    PRINT_DIMSET(dimViscosity);
    PRINT_DIMSET(dimDynamicViscosity);


    Info<< nl << "Operators" << nl << nl;

    // Operators
    {
        Info<< "dimLength/dimTime  " << (dimLength/dimTime) << nl;
        Info<< "1/dimTime  " << (dimless/dimTime) << nl;
        Info<< "-dimTime  " << (-dimTime) << nl;
        Info<< "~dimTime  " << (~dimTime) << nl;

        Info<< "sqr(dimVelocity) " << sqr(dimVelocity) << nl;
        Info<< "pow2(dimVelocity) " << pow2(dimVelocity) << nl;
        Info<< "sqrt(dimVelocity) " << sqrt(dimVelocity) << nl;
        Info<< "pow025(dimVelocity) " << pow025(dimVelocity) << nl;
        Info<< "sqrt(sqrt(dimVelocity)) " << sqrt(sqrt(dimVelocity)) << nl;
    }


    unsigned nFail = 0;

    nFail += checkDimensions
    ({
        { true, dimless, dimless },
        { false, dimless, dimPressure }
    });


    if (nFail)
    {
        Info<< nl << "failed " << nFail << " tests" << nl;
        return 1;
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
