/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    Test-Function1

Description
    Tests Function1

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "Function1.H"
#include "scalarIndList.H"
#include "scalarField.H"
#include "IOdictionary.H"
#include "linearInterpolationWeights.H"
#include "splineInterpolationWeights.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    const word dictName("function1Properties");

    argList::addBoolOption("all", "Test all functions in function1Properties");

    argList::addArgument("function1");
    argList::addArgument("...");
    argList::addArgument("functionN");
    argList::noMandatoryArgs();

    #include "setRootCase.H"
    #include "createTime.H"

    {
        scalarField samples({0, 1, 2, 3});

        scalarField values(4, scalar(1));

        linearInterpolationWeights interpolator
        //splineInterpolationWeights interpolator
        (
            samples
        );
        labelList indices;
        scalarField weights;

        interpolator.integrationWeights(1.1, 1.2, indices, weights);
        Pout<< "indices:" << indices << nl
            << "weights:" << weights << nl;

        scalar baseSum = interpolator.weightedSum
        (
            weights,
            scalarUIndList(values, indices)
        );
        Pout<< "baseSum=" << baseSum << nl << nl << endl;

        // interpolator.integrationWeights(-0.01, 0, indices, weights);
        // scalar partialSum = interpolator.weightedSum
        // (
        //     weights,
        //     scalarUIndList(values, indices)
        // );
        // Pout<< "partialSum=" << partialSum << nl << nl << endl;

        // interpolator.integrationWeights(-0.01, 1, indices, weights);
        // //Pout<< "samples:" << samples << endl;
        // //Pout<< "indices:" << indices << endl;
        // //Pout<< "weights:" << weights << endl;
        // scalar sum = interpolator.weightedSum
        // (
        //     weights,
        //     scalarUIndList(values, indices)
        // );
        // Pout<< "integrand=" << sum << nl << nl << endl;
    }

    if (args.found("all") || args.size() > 1)
    {
        #include "setConstantRunTimeDictionaryIO.H"

        #if (OPENFOAM > 2212)
        dictionary propsDict(IOdictionary::readContents(dictIO));
        #else
        dictionary propsDict(static_cast<dictionary&&>(IOdictionary(dictIO)));
        #endif

        const scalarField xvals(propsDict.lookup("x"));

        Info<< "Entries" << flatOutput(propsDict.toc()) << nl << nl;

        Info<< "Inputs" << nl
            << "    x = " << xvals << nl
            << endl;

        DynamicList<word> functionNames;

        auto nameFilter = [](const word& val)
        {
            return !(val == "x" || val.ends_with("Coeffs"));
        };

        if (args.found("all"))
        {
            for (const word& f : propsDict.toc())
            {
                if (nameFilter(f))
                {
                    functionNames.push_back(f);
                }
            }
        }
        else
        {
            for (label argi=1; argi < args.size(); ++argi)
            {
                functionNames.push_back(args[argi]);
            }
        }

        for (const word& funName : functionNames)
        {
            auto function1 = Function1<scalar>::New(funName, propsDict);

            // Info<< "Data entry type: " << function1().type() << nl;
            Info<< "////" << nl;
            function1().writeData(Info);
            Info<< nl;

            Info<< "Values" << nl;
            for (const scalar& x : xvals)
            {
                Info<< "    f(" << x << ") = " << function1().value(x) << nl;
            }

            Info<< endl;
        }
    }

    return 0;
}


// ************************************************************************* //
