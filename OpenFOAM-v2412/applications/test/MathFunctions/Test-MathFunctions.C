/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    Test-MathFunctions

Description
    Tests for \c Math namespace member functions
    using \c doubleScalar base type.

\*---------------------------------------------------------------------------*/

#include "MathFunctions.H"
#include "mathematicalConstants.H"
#include "IOmanip.H"
#include "TestTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Execute each member functions of Foam::Math::, and print output
template<class Type>
void test_member_funcs(Type)
{
    Random rndGen(1234);
    const label numberOfTests = 10000;
    const scalar tolerance = 1e-3;

    for (label i = 0; i < numberOfTests; ++i)
    {
        Info<< "# Inverse error functions:" << endl;
        const Type x1 = rndGen.sample01<Type>();
        const Type x2 = rndGen.sample01<Type>();
        Info<< "    # Operands:" << nl
            << "      x1 = " << x1 << nl
            << "      x2 = " << x2 << endl;

        cmp
        (
            "       # erf(erfinv(x1)) = x1 = ",
            x1,
            Foam::erf(Foam::Math::erfInv(x1)),
            tolerance
        );

        cmp
        (
            "       # erf(erfinv(-x2)) = -x2 = ",
            x2,
            Foam::erf(Foam::Math::erfInv(x2)),
            tolerance
        );

        Info<< "# Incomplete gamma functions:" << endl;
        const Type a = rndGen.position(SMALL, Type(100));
        const Type x = rndGen.position(Type(0), Type(100));
        Info<< "    # Operands:" << nl
            << "      a = " << a << nl
            << "      x = " << x << endl;

        cmp
        (
            "       # incGammaRatio_Q(a,x) + incGammaRatio_P(a,x) = 1 = ",
            Foam::Math::incGammaRatio_Q(a,x) + Foam::Math::incGammaRatio_P(a,x),
            scalar(1),
            tolerance
        );

        cmp
        (
            "       # incGamma_Q(1, x) = exp(-x) = ",
            Foam::Math::incGamma_Q(1, x),
            Foam::exp(-x),
            tolerance
        );

        cmp
        (
            "       # incGamma_Q(0.5, x) = sqrt(pi)*erfc(sqrt(x)) = ",
            Foam::Math::incGamma_Q(0.5, x),
            Foam::sqrt(Foam::constant::mathematical::pi)
           *Foam::erfc(Foam::sqrt(x)),
            tolerance
        );

        cmp
        (
            "       # incGamma_P(1,x) = 1 - exp(x) = ",
            Foam::Math::incGamma_P(1, x),
            1 - Foam::exp(-x),
            tolerance
        );

        cmp
        (
            "       # incGamma_P(0.5, x) = sqrt(pi)*erf(sqrt(x)) = ",
            Foam::Math::incGamma_P(0.5, x),
            Foam::sqrt(Foam::constant::mathematical::pi)
           *Foam::erf(Foam::sqrt(x)),
            tolerance
        );
    }
}


// Do compile-time recursion over the given types
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
run_tests(const std::tuple<Tp...>& types, const List<word>& typeID){}


template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
run_tests(const std::tuple<Tp...>& types, const List<word>& typeID)
{
    Info<< nl << "    ## Test member functions: "<< typeID[I] <<" ##" << nl;
    test_member_funcs(std::get<I>(types));

    run_tests<I + 1, Tp...>(types, typeID);
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main()
{
    Info<< setprecision(15);

    const std::tuple<doubleScalar> types
    (
        std::make_tuple(Zero)
    );

    const List<word> typeID
    ({
        "doubleScalar"
    });

    run_tests(types, typeID);


    if (nFail_)
    {
        Info<< nl << "        #### "
            << "Failed in " << nFail_ << " tests "
            << "out of total " << nTest_ << " tests "
            << "####\n" << endl;
        return 1;
    }

    Info<< nl << "        #### Passed all " << nTest_ <<" tests ####\n" << endl;
    return 0;
}


// ************************************************************************* //
