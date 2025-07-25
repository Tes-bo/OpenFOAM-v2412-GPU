/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Test functionality of DirLister

\*---------------------------------------------------------------------------*/

#include "DirLister.H"
#include "fileNameList.H"
#include "wordRes.H"
#include "predicates.H"
#include "FlatOutput.H"
#include "error.H"
#include "stringOps.H"
#include "scalar.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    {
        Info<< nl
            << "List items" << nl
            << "~~~~~~~~~~" << nl;
        for (const word& item : DirLister("."))
        {
            Info<< "    " << item << nl;
        }
    }

    {
        Info<< nl
            << "List files" << nl
            << "~~~~~~~~~~" << nl;
        for (const word& item : DirLister::files("."))
        {
            Info<< "    " << item << nl;
        }
    }

    {
        Info<< nl
            << "List dirs" << nl
            << "~~~~~~~~~~" << nl;
        for (const word& item : DirLister::dirs("."))
        {
            Info<< "    " << item << nl;
        }
    }

    {
        Info<< nl
            << "List files - filtered" << nl
            << "~~~~~~~~~~" << nl;
        for
        (
            const word& item
          : DirLister::files(".").where
            (
                [](const word& val){ return val.starts_with('T'); }
            )
        )
        {
            Info<< "    " << item << nl;
        }
    }

    {
        Info<< nl
            << "List dirs - filtered" << nl
            << "~~~~~~~~~~" << nl;

        for
        (
            const word& item
          : DirLister::dirs(".").where(regExp("Ma.*"))
        )
        {
            Info<< "    " << item << nl;
        }
    }

    {
        Info<< nl
            << "List dirs - filtered" << nl
            << "~~~~~~~~~~" << nl;

        for
        (
            const word& item
          : DirLister::dirs(".").where(predicates::always())
        )
        {
            Info<< "    " << item << nl;
        }
    }


    {
        Info<< nl
            << "List items" << nl
            << "~~~~~~~~~~" << nl
            << DirLister(".").list<fileName>() << nl;
    }

    {
        Info<< nl
            << "List files - filtered" << nl
            << "~~~~~~~~~~" << nl
            <<  DirLister(".").list<fileName>
                (
                    [](const word& val){ return val.starts_with('D'); },
                    false
                )
            << nl;
    }

    {
        Info<< nl
            << "List files - filtered" << nl
            << "~~~~~~~~~~" << nl;

        wordRes relist
        ({
            wordRe("processors"),
            wordRe("processor[0-9][0-9]*", wordRe::REGEX)
        });

        Info<<"matcher: " << flatOutput(relist) << endl;

        for (const word& item : DirLister::dirs(".").where(relist))
        {
            Info<< "=>    " << item << nl;
        }

        Info<< "dirList: "
            << flatOutput
               (
                   DirLister::dirs(".").csorted<fileName>(relist)
               ) << nl;
    }


    {
        Info<< nl
            << "List time dirs" << nl
            << "~~~~~~~~~~" << nl;

        for
        (
            const word& item
          : DirLister::dirs(".").where
            (
                [](const word& val)
                {
                    scalar s;
                    return readScalar(val, s) || val == "constant";
                }
            )
        )
        {
            Info<< "=>    " << item << nl;
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
