/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Test minMax

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Time.H"
#include "BitOps.H"
#include "HashOps.H"
#include "ListOps.H"
#include "scalarField.H"
#include "complexField.H"
#include "MinMax.H"
#include "dimensionedScalar.H"

using namespace Foam;


template<class T>
Ostream& printInfo(const MinMax<T>& range)
{
    Info<< range << " good=" << range.good() << " span=" << range.span();

    return Info;
}


template<class T>
void testUniformField(const T& val)
{
    constexpr label N = 10;

    //    Field<T> fld(N, val);
    List<T> fld(N, val);

    Info<< "field:   " << fld << nl
        << "min/max: " << minMaxMag(fld) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    Info<< "Test min/max " << nl;

    Info<<"Construct null: ";
    printInfo(MinMax<scalar>()) << nl;

    Info<<"Construct single value : ";
    printInfo(MinMax<scalar>(15)) << nl;

    Info<<"Construct zero : ";
    printInfo(MinMax<scalar>(Zero)) << nl;

    Info<<"Construct zero_one : ";
    printInfo(MinMax<scalar>(zero_one{})) << nl;

    Info<<"Construct range : ";
    printInfo(MinMax<scalar>(1, 20)) << nl;

    Info<<"A 0-1 scalar range : ";
    printInfo(scalarMinMax::zero_one()) << nl;

    Info<<"A 0-1 vector range : ";
    printInfo(MinMax<vector>::zero_one()) << nl;

    {
        vector a(0, 1, 20);
        vector b(2, 1, 0);
        vector c(4, 10, 12);

        Info<< "vectors:"
            << " a = " << a
            << " b = " << b
            << " c = " << c << nl;

        Info<< "min max    = " << min(max(a, b), c) << nl;
        Info<< "range clamp= " << MinMax<vector>(b, c).clamp(a) << nl;
        Info<< "clamp      = " << clamp(a, b, c) << nl;

        Info<< "clamp 0/1  = " << clamp(a, vector::zero, vector::one) << nl;
    }

    // Scalar promotion
    Info<< "clamp (scalar) = " << clamp(15.0, -1, 1) << nl;


    {
        scalarMinMax range1(10, 20);
        scalarMinMax range2(40, 50);
        Info<< range1 << " + " << range2 << " = " << (range1 + range2) <<nl;
    }

    {
        Info<<"Dimensioned range : "
            << dimensioned<scalarMinMax>("velrange", dimVelocity, {1, 20})
            << nl;

        dimensioned<scalarMinMax> range1("a", dimVelocity, {10, 20});
        dimensioned<scalarMinMax> range2("b", dimVelocity, {40, 50});

        Info<<"Dimensioned range : " << (range1 + range2) << endl;
    }

    Info<<"Centre value for (10 250) = "
        << scalarMinMax(10, 250).centre() << nl;


    {
        Info<<"compare - similar definition as per std::string compare" << nl;

        string str1("abc");
        Info<< "For string=" << str1 << nl
            << "    compare(\"aaa\") = " << str1.compare("aaa") << nl
            << "    compare(\"abz\") = " << str1.compare("abz") << nl;


        MinMax<scalar> range(10, 20);

        Info<< "For range=" << range << nl
            << "    compare(5) = " << range.compare(5) << nl
            << "    compare(15) = " << range.compare(15) << nl
            << "    compare(25) = " << range.compare(25) << nl;


        Info<< "Binary comparisons" << nl;

        Info<< "(5 < range) = " << (5 < range) << nl
            << "(25 > range) = " << (25 >= range) << nl
            << "(12 <= range) = " << (12 <= range) << nl
            << "(12 >= range) = " << (12 >= range) << nl;
    }


    scalarField values1
    (
        List<scalar>({3, 10, -11, 85, 300})
    );

    values1 *= (Pstream::myProcNo()+1);

    Pout<<"min-max of " << flatOutput(values1) << " = "
        << minMax(values1) << endl;

    // Construct from values
    MinMax<scalar> minmax1(values1);
    printInfo(minmax1) << nl;

    // Reset and add values
    minmax1.clear();

    minmax1 += values1;
    Pout<<"range: " << minmax1 << endl;


    Info<< "Reduced: "<< returnReduce(minmax1, plusOp<scalarMinMax>()) << nl;
    Info<< "Reduced: "<< returnReduce(minmax1, minMaxOp<scalar>()) << nl;

    // Info<< "gMinMax: "<< gMinMax(values1v) << nl;

    vectorField values1v
    (
        ListOps::create<vector>
        (
            values1,
            [](const scalar s) { return vector(s, 2*s, -2*s); }
        )
    );

    Info<< "gMinMax: " << gMinMax(values1v) << nl;
    Info<< "gMinMaxMag: " << gMinMaxMag(values1v) << nl;

    {
        MinMax<scalar> limiter(10, 200);
        clampOp<scalar> clipper(limiter);
        // clampOp<scalar> clipper(zero_one{});
        // clampOp<scalar> clipper(10, 200);

        Info<< nl
            << "Test clipping limiter: " << limiter << nl
            << "values : " << flatOutput(values1) << nl;

        Info<< "Subset mask: "
            << ListOps::create<bool, scalar>(values1, limiter)
            << nl;

        Info<< "Subset = " << subsetList(values1, limiter) << nl;
        Info<< "Subset = "
            << IndirectList<scalar>::subset_if(values1, limiter) << nl;

        Info<< nl << "test clip() with limiter: " << limiter << nl;
        for (const scalar& val : values1)
        {
            Info<< "clipped : " << val << " = "
                << clip(val, limiter)
                << " or " << clip(val, limiter)
                << " or " << clipper(val) << nl;
        }

        Info<< nl << "test clip(Field) with limiter: " << limiter << nl;
        Info<< "clipped : " << clip(values1, limiter) << nl;

        scalarField values2(values1);

        Info<< nl << "inplace clip" << nl;

        Info<< "before " << flatOutput(values2) << nl;

        // Too much clutter
        // for (scalar& val : values2)
        // {
        //     clampEqOp<scalar>{limiter}(val);
        // }

        for (scalar& val : values2)
        {
            val = clipper(val);
        }

        Info<< "after: " << flatOutput(values2) << nl;

        Info<< nl << "For list: " << flatOutput(values1) << nl
            << " minMax    : " << minMax(values1) << nl
            << " minMaxMag : " << minMaxMag(values1)
            << " = " << mag(minMaxMag(vector(1, 2, 3)))
            << nl;
    }


    // Hashed values and reduction
    {
        HashTable<scalarMinMax> hashed;

        hashed.insert("zero", scalarMinMax(Zero));
        hashed("null");
        hashed("U") += values1;


        Pout<< "hashed: " << hashed << nl;

        Pstream::mapCombineReduce(hashed, plusEqOp<scalarMinMax>());

        Info<< "reduced: " << hashed << nl;


        // Prune invalid
        hashed.filterValues(emptyOp<scalarMinMax>(), true);

        Info<< "filtered: " << hashed << nl;
    }


    // Min/max of uniform fields
    {
        testUniformField<scalar>(100);
        // testUniformField<complex>(complex(100, 0));
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
