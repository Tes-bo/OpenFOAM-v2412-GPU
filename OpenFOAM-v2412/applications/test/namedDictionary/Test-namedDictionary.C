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
    Test-namedDictionary

Description
    Test handling of keyType/dictionary

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "namedDictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::addArgument("file1 .. fileN");
    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        InfoErr<< "Provide a file or files to test" << nl;
    }
    else
    {
        for (label argi=1; argi < args.size(); ++argi)
        {
            const auto dictFile = args.get<fileName>(argi);
            IFstream ifs(dictFile);

            dictionary dict(ifs);

            IOobject::writeDivider(Info) << nl;

            for (const entry& dEntry : dict)
            {
                if (dEntry.isStream())
                {
                    List<namedDictionary> list(dEntry.stream());

                    Info<< "input: " << dEntry << nl
                        << "list: " << list << nl;
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
