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

Application
    Test-vector

Description
    Some simple tests for vector

\*---------------------------------------------------------------------------*/

#include "vectorField.H"
#include "IOstreams.H"
#include "Random.H"
#include <algorithm>
#include <random>

using namespace Foam;

void printInfo(const vector& vec)
{
    Info<< "vector : " << vec << nl
        << "magSqr : " << magSqr(vec) << nl;

    Info<< "component"
        << " max:" << cmptMax(vec)
        << " sum:" << cmptSum(vec)
        << " prod:" << cmptProduct(vec)
        << " mag:" << cmptMag(vec)
        << " magSqr:" << cmptMagSqr(vec)
        << nl << nl;
}


void doTest(vector& vec1, vector& vec2)
{
    Info<<"normalised(vector1): " << normalised(vec1) << nl;
    Info<<"vector1: " << vec1 << nl;
    vec1.normalise();
    Info<<"normalised: " << vec1 << nl;

    vector vecsmall((1e-100 * vec1));

    Info<<"small: " << vecsmall << nl;
    Info<<"small normalised: " << normalised(vecsmall) << nl;
    Info<<"small diff: " << (vecsmall.normalise() - vec1) << nl;

    vec1 *= 4.0;

    Info<< "scalar mult: " << vec1 << nl;
    Info<< "addition   : " << (vec1 + vec1) << nl;

    printInfo(vec1);
    printInfo(vec2);

    Info<< "vector: " << vec1 << nl
        << "vector: " << vec2 << nl
        << "   min: " << min(vec1, vec2) << nl
        << "  dist: " << vec1.dist(vec2) << ' ' << mag(vec1 - vec2) << nl
        << "dist^2: " << vec1.distSqr(vec2) << ' ' << magSqr(vec1 - vec2) << nl
        << nl;
}


template<class VecSpace>
void testIterator(const VecSpace& vs)
{
    Info<< "size: " << vs.size() << " for:";
    for (const auto& val : vs)
    {
        Info<< " " << val;
    }
    Info<< nl;
}


template<class VecSpace>
void testData(const VecSpace& vs)
{
    Info<< "size: " << vs.size() << " for:";

    const auto* data = vs.cdata();
    for
    (
        const auto* endData = data + VecSpace::nComponents;
        data != endData;
        ++data
    )
    {
        Info<< " " << *data;
    }
    Info<< nl;
}


template<class Type>
void testNormalise(Field<Type>& fld)
{
    Info<< nl << pTraits<Type>::typeName << " Field" << nl
        << "  orig: " << fld << nl;
    fld.normalise();
    Info<< "  norm: " << fld << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<<"normalised: " << vector(1,2,3).normalise() << nl;
    Info<<"normalised: " << vector::uniform(VSMALL).normalise() << nl;
    Info<<"normalised: " << vector::uniform(ROOTVSMALL).normalise() << nl;

    {
        vector vec1(0.5, 0.5, 0.5);
        vector vec2(0.5, 0.51, -0.5);

        doTest(vec1, vec2);

        testIterator(vec1);
        testIterator(vec2);
        testData(vec2);

        // Use STL algorithm(s)

        std::sort(vec2.begin(), vec2.end());
        Info<< "sorted: " << vec2 << nl;

        std::shuffle(vec2.begin(), vec2.end(), std::default_random_engine());
        Info<< "shuffled: " << vec2 << nl;

        // Vectors with some identical components
        vectorField vectors
        ({
            {1.1, 2.2, 3.3 },
            {2.2, 3.3, 4.4 },
            {-1.1, 2.2, 3.3 },
            {-2.2, 3.3, 4.4 },

            {-1.1, -2.2, 3.3 },
            {-2.2, -3.3, 4.4 },

            {-1.1, -2.2, -3.3 },
            {-2.2, -3.3, -4.4 },
            {-3.3, 2.1, 12 },
            {3.145, 1.6, 2 },

            {0, 0, 0}
        });

        shuffle(vectors);

        Info<< "initial vectors: ";
        vectors.writeList(Info, 1) << nl;

        Foam::sort(vectors);
        Info<< "regular sort:";
        vectors.writeList(Info, 1) << nl;

        std::sort(vectors.begin(), vectors.end(), vector::less_xyz);
        Info<< "sorted xyz:";
        vectors.writeList(Info, 1) << nl;

        std::sort(vectors.begin(), vectors.end(), vector::less_yzx);
        Info<< "sorted yzx:";
        vectors.writeList(Info, 1) << nl;

        std::sort(vectors.begin(), vectors.end(), vector::less_zxy);
        Info<< "sorted zxy:";
        vectors.writeList(Info, 1) << nl;

        vectorField vecCmptMag1(cmptMag(vectors));
        Info<< "cmptMag:";
        vecCmptMag1.writeList(Info, 1) << nl;

        vectorField vecCmptMag2(vectors.size());
        vecCmptMag2.replace(vector::X, mag(vectors.component(vector::X)));
        vecCmptMag2.replace(vector::Y, mag(vectors.component(vector::Y)));
        vecCmptMag2.replace(vector::Z, mag(vectors.component(vector::Z)));

        Info<< "cmptMag:";
        vecCmptMag2.writeList(Info, 1) << nl;

        Info<< "mult:";
        cmptMultiply(vecCmptMag2, vecCmptMag2, vector(2,3,4));
        vecCmptMag2.writeList(Info, 1) << nl;
    }
    // Basic tests for fields
    {
        scalarField sfld1
        ({
            0.0,
            1.0,
            2.0
        });

        testNormalise(sfld1);

        Field<floatVector> vfld1
        ({
            floatVector(0.0, 1.0, 2.0),
            floatVector(0, 0, floatScalarSMALL),
            floatVector(0, 0, floatScalarVSMALL),
            floatVector(0, 2, 1),
        });

        testNormalise(vfld1);

        Field<doubleVector> vfld2
        ({
            doubleVector(0.0, 1.0, 2.0),
            doubleVector(0, 0, doubleScalarSMALL),
            doubleVector(0, 0, doubleScalarVSMALL),
            doubleVector(0, 2, 1),
        });

        testNormalise(vfld2);
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
