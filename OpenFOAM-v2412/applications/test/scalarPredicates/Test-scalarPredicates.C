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

Application
    Test-scalarPredicates

Description
    Simple tests using predicates for scalars

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "labelList.H"
#include "scalarList.H"
#include "scalarPredicates.H"
#include "FlatOutput.H"
#include "Tuple2.H"
#include "SpanStream.H"
#include "ops.H"
#include "bitSet.H"

using namespace Foam;


void doTest(const scalarList& values, const predicates::scalars& accept)
{
    // Also tests that output is suppressed
    Info<<"Have: " << accept.size() << " predicates" << accept << endl;
    Info<<"values: " << flatOutput(values) << endl;

    for (const scalar& value : values)
    {
        if (accept.match(value))
        {
            Info<< "matched: " << value << " at "
                << flatOutput(accept.matching(value)) << nl;
        }
    }

    labelList matches = accept.matching(values);
    Info<< "values matched at positions: " << flatOutput(matches) << nl;
}


template<class Predicate>
void testPredicate(const scalarList& values, const Predicate& pred)
{
    bitSet matches;

    label i=0;

    for (const scalar& value : values)
    {
        if (pred(value))
        {
            matches.set(i);
        }

        ++i;
    }

    IndirectList<scalar> matched(values, matches.toc());

    Info<< "matched: " << flatOutput(matched.addressing())
        << " = " << flatOutput(matched) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    scalarList values
    ({
        -10,
        10,
        0,
        3.145,
        1000.56,
        1e5,
    });


    Info<< nl << "From a mixed list of entries" << nl;
    {
        predicates::scalars
            accept
            ({
                predicates::scalars::lessOp(10),
                predicates::scalars::greaterOp(100),
                predicates::scalars::orOp
                (
                    predicates::scalars::lessOp(-5),
                    predicates::scalars::greaterOp(100)
                ),
                predicates::scalars::orOp
                (
                    [](const scalar& val){ return val < 5; },
                    predicates::scalars::greaterOp(100)
                ),

                [](const scalar& val){ return val < -8; },

                // Rather wordy, word normally not be called manually
                predicates::scalars::operation("le", -9),
            });

        doTest(values, accept);
    }


    Info<< nl << "Construct from list input" << nl;
    {
        List<Tuple2<word, scalar>> entries
        ({
            {"less", 10},
            {"greater", 100},
            // Not possible >  ((less -5) or (greater 100))
            {"less", -8},
            {"le", -9},
        });

        predicates::scalars accept(entries);

        doTest(values, accept);
    }

    Info<< nl << "Construct from initializer_list input" << nl;
    {
        predicates::scalars accept
        ({
            {"less", 10},
            {"greater", 100},
            // Not possible >  ((less -5) or (greater 100))
            {"less", -8},
            {"le", -9},
        });

        doTest(values, accept);
    }

    Info<< nl << "Construct from Istream" << nl;
    {
        ICharStream is("((less 10) (greater 100) (less -8) (le -9))");
        predicates::scalars accept(is);

        doTest(values, accept);

        // change some location
        accept[0] = predicates::scalars::greaterOp(1000);

        Info<< nl << "Reset some values" << nl;
        doTest(values, accept);
    }


    Info<< nl << "Test with ops" << nl;
    Info<<"values: " << flatOutput(values) << endl;
    {
        testPredicate(values, lessOp1<scalar>(10));
        testPredicate(values, greaterOp1<scalar>(100));

        // Also with dissimilar type
        testPredicate(values, lessEqualOp1<label>(0));
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
