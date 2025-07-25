/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Test-codeStream

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

    Info<< nl
        << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
        << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
        << endl;

    if (args.size() <= 1)
    {
        Info<<"specify dictionaries to test\n";
    }
    else
    {
        IOobject::writeDivider(Info);
        for (label argi=1; argi < args.size(); ++argi)
        {
            const auto dictFile = args.get<fileName>(argi);
            IFstream is(dictFile);

            dictionary dict(is);

            Info<< dict << endl;
            for (const entry& e : dict)
            {
                if (e.isStream())
                {
                    Info<< e.keyword() << " :: "
                        << e.stream().toString() << nl;
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
