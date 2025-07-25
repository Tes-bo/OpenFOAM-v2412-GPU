/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Simple field tests

    Test use of Kahan/Neumaier to extend precision for when running SPDP
    mode. Conclusion is that it is easier/quicker to run these summation
    loops as double precision (i.e. solveScalar).

\*---------------------------------------------------------------------------*/

#include "primitiveFields.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp, class ResultType>
void sumNeumaier
(
    const UList<Type>& vals,
    const CombineOp& cop,
    ResultType& result
)
{
    // Neumaier version of Kahan
    ResultType sum = Zero;
    ResultType c = Zero;
    for (const Type& vali : vals)
    {
        ResultType val;
        cop(val, vali);

        const ResultType t = sum + val;

        for
        (
            direction cmpt = 0;
            cmpt < pTraits<ResultType>::nComponents;
            cmpt++
        )
        {
            if (mag(sum[cmpt]) >= mag(val[cmpt]))
            {
                // If sum is bigger, low-order digits of input[i] are lost.
                c[cmpt] += (sum[cmpt] - t[cmpt]) + val[cmpt];
            }
            else
            {
                // Else low-order digits of sum are lost.
                c[cmpt] += (val[cmpt] - t[cmpt]) + sum[cmpt];
            }
        }
        sum = t;
    }
    result = sum + c;
}


template<class CombineOp, class ResultType>
void sumNeumaier
(
    const UList<scalar>& vals,
    const CombineOp& cop,
    ResultType& result
)
{
    // Neumaier version of Kahan
    ResultType sum = Zero;
    ResultType c = Zero;
    for (const scalar vali : vals)
    {
        ResultType val;
        cop(val, vali);

        const ResultType t = sum + val;
        if (mag(sum) >= mag(val))
        {
            // If sum is bigger, low-order digits of input[i] are lost.
            c += (sum - t) + val;
        }
        else
        {
            // Else low-order digits of sum are lost.
            c += (val - t) + sum;
        }
        sum = t;
    }
    result = sum + c;
}


template<class Type>
Type mySum(const UList<Type>& f)
{
    typedef typename Foam::typeOfSolve<Type>::type solveType;

    solveType Sum = Zero;

    if (f.size())
    {
        sumNeumaier(f, eqOp<solveType>(), Sum);
    }

    return Type(Sum);
}


//- The sumSqr always adds only positive numbers. Here there can never be any
//  cancellation of truncation errors.
template<class Type>
typename outerProduct1<Type>::type mySumSqr(const UList<Type>& f)
{
    typedef typename outerProduct1<solveScalar>::type prodType;

    prodType result = Zero;

    if (f.size())
    {
        sumNeumaier(f, eqSqrOp<prodType>(), result);
    }

    return result;
}


// Main program:

int main(int argc, char *argv[])
{
    scalarField sfield(10, one{});

    forAll(sfield, i)
    {
        sfield[i] = (i % 4) ? i : 0;
    }

    Info<< "scalarField: " << sfield << nl;
    sfield.negate();

    Info<< "negated: " << sfield << nl;

    // Does not compile (ambiguous)
    // boolField lfield(10, one{});

    boolField lfield(10, true);

    forAll(lfield, i)
    {
        lfield[i] = (i % 4) ? i : 0;
    }

    Info<< "boolField: " << lfield << nl;
    lfield.negate();

    Info<< "negated: " << lfield << nl;


    // Summation (compile in SPDP)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Pout.precision(16);
    Sout.precision(16);

    const scalar SMALLS(1e-6);
    const scalar GREATS(1e6);

    // scalarField summation
    {
        scalarField sfield(10, SMALLS);

        sfield[8] = GREATS;
        sfield[9] = -sfield[8];
        Info<< "scalarField:" << sfield.size() << nl
            << "    sum      :" << sum(sfield) << nl
            << "    corrected:" << mySum(sfield) << endl;
    }
    // vectorField summation
    {
        vectorField vfield(10, vector::uniform(SMALLS));

        vfield[8] = vector::uniform(GREATS);
        vfield[9] = -vfield[8];
        Info<< "vectorField:" << vfield.size() << nl
            << "    sum      :" << sum(vfield) << nl
            << "    corrected:" << mySum(vfield) << endl;
    }
    // sphericalTensorField summation
    {
        sphericalTensorField tfield(10, sphericalTensor(SMALLS));

        tfield[8] = sphericalTensor(GREATS);
        tfield[9] = -tfield[8];
        Info<< "sphericalTensorField:" << tfield.size() << nl
            << "    sum      :" << sum(tfield) << nl
            << "    corrected:" << mySum(tfield) << endl;
    }
    // symmTensorField summation
    {
        symmTensorField tfield(10, SMALLS*symmTensor::I);

        tfield[8] = GREATS*symmTensor::I;
        tfield[9] = -tfield[8];
        Info<< "symmTensorField:" << tfield.size() << nl
            << "    sum      :" << sum(tfield) << nl
            << "    corrected:" << mySum(tfield) << endl;
    }
    // tensorField summation
    {
        tensorField tfield(10, SMALLS*tensor::I);

        tfield[8] = GREATS*tensor::I;
        tfield[9] = -tfield[8];
        Info<< "tensorField:" << tfield.size() << nl
            << "    sum      :" << sum(tfield) << nl
            << "    corrected:" << mySum(tfield) << endl;
    }


    return 0;
}


// ************************************************************************* //
