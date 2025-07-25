/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Test-charList

Description
    Some test of UList, List for characters

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "IOstreams.H"
#include "messageStream.H"

#include "charList.H"
#include "labelList.H"
#include "SpanStream.H"
#include "ListOps.H"
#include "SubList.H"
#include "FlatOutput.H"

#include <numeric>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"

    // Info<< "Known compound tokens: "
    //     << token::compound::emptyConstructorTablePtr_->sortedToc() << nl;

    OCharStream ostr;

    {
        List<char> alphabet(64);
        std::iota(alphabet.begin(), alphabet.end(), '@');

        alphabet.writeEntry("alphabet", Info);

        ostr << alphabet;

        Info<< "wrote: " << ostr.str() << nl;
    }

    {
        // ICharStream istr(ostr.release());
        ISpanStream istr(ostr.view());
        List<char> alphabet(istr);

        Info<< "re-read: " << alphabet << nl;

        // Can assign zero?
//Fails:        alphabet = char(Zero);
        alphabet = Foam::zero{};

        // alphabet = '@';
        Info<< "blanked: " << alphabet << nl;
    }

    Info<< "\nEnd\n" << nl;

    return 0;
}

// ************************************************************************* //
