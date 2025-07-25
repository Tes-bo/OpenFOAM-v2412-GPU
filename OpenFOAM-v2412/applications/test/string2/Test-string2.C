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
    Test some string functionality

\*---------------------------------------------------------------------------*/

#include "fieldExprFwd.H"
#include "int.H"
#include "uint.H"
#include "scalar.H"
#include "Switch.H"
#include "fileName.H"
#include "stringList.H"
#include "stringOps.H"
#include "dictionary.H"
#include "IOstreams.H"
#include "OSspecific.H"


using namespace Foam;

void testCommentStripping(const std::string& s)
{
    Info<< "input" << nl
        << "========" << nl
        << s << nl
        << "========" << nl;

    Info<< "output" << nl
        << "========" << nl
        << stringOps::removeComments(s) << nl
        << "========" << nl << nl;
}


void testNumericEvaluation(const std::string& s)
{
    Info<< "input" << nl
        << "========" << nl
        << s << nl
        << "========" << nl;

    const bool oldThrowingIOErr = FatalIOError.throwing(true);

    try
    {
        std::string expanded(stringOps::expand(s));

        if (expanded == s)
        {
            Info<< "DID NOT EXPAND" << nl;
        }
        else
        {
            Info<< "output" << nl
                << "========" << nl
                << expanded << nl;
        }
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "Expand triggered FatalIOError:"
            << err.message().c_str() << nl;
    }

    Info<< "========" << nl << nl;

    FatalIOError.throwing(oldThrowingIOErr);
}


// Test findTrim - uses '<' and '>' as pseudo-placeholders
void testFindTrim(const std::string& s)
{
    auto pos = s.find('<');
    auto len = s.find('>');

    // Conform with expected value for substr

    if (pos == std::string::npos) pos = 0; else ++pos;
    if (len != std::string::npos) len = (len - pos);

    const auto pts = stringOps::findTrim(s, pos, len);

    Info<< "input" << nl
        << "========" << nl
        << s << nl
        << "pos=" << pos << " len=" << int(len) << nl
        << s.substr(pos, len) << nl
        << "trim=" << pts.first << " to " << pts.second << nl
        << s.substr(pts.first, pts.second-pts.first) << nl
        << "========" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Test comment stripping
    {
        Info<< nl << "Test comment stripping" << nl;

        for
        (
            const auto& cstr
          :
            {
                "/String without comments/",
                "Removed some/* C-comments */ / C comments",
                "Removed some//C++ comments\n / C++ comments",
                "Partly degenerate C comment </*/ C-comment...",
                "Truncated C comment </* C-comment...",
                "Truncated C++ comment <// C++ comment...",
            }
        )
        {
            testCommentStripping(cstr);
        }
    }

    // Test numeric
    {
        Info<< nl << "Test numeric evaluation" << nl;

        // expressions::fieldExpr::debug = 2;

        for
        (
            const auto& cstr
          :
            {
                "My value <${{ round(100 / 15) }}> as int",
                "sqrt(2) = (${{ sqrt(2) }})",
                "sqrt(2) = (${{ sqrt(2) }/* Truncated */",
                "sqrt(2) = (${{ sqrt(2) * foo() }})/* bad expr */",
                "huge = (${{ sqrt(123E+5000) }})/* range error */",

                "vector=${{   5 * vector(1,2,3)   }}=",
                "empty=${{    }}=",

                "vector=${{   5 * vector(1, 2, 3) ^ vector(4, 5, 6) }}=",
// NOT YET WORKING: "vector=${{   5 * [1 2 3 ] ^ [ 4 5 6 ] }}=",
            }
        )
        {
            testNumericEvaluation(cstr);
        }
    }


    // Test findTrim - uses '<' and '>' as pseudo-placeholders
    {
        Info<< nl << "Test findTrim" << nl;

        for
        (
            const auto& cstr
          :
            {
                "",     // Empty
                "   \n\t   ",   // whitespace only
                "   Leading and trailing   ",
                "   Alt end and trai>ling   ",
                "   Alt start<   space      ",
            }
        )
        {
            testFindTrim(cstr);
        }
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
