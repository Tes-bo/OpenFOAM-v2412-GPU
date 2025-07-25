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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-dictionaryCopy

Description
    Test copying a dictionary with filtering

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "dictionaryContent.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

    const wordRes allow;
    const wordRes deny
    ({
        wordRe("boundary.*", wordRe::REGEX)
    });

    Info<< nl
        << "allow: " << flatOutput(allow) << nl
        << "deny:  " << flatOutput(deny) << nl << nl;

    if (args.size() <= 1)
    {
        const string dictFile = "testDictCopy";
        IFstream is(dictFile);

        dictionary input(is);

        dictionary copied
        (
            dictionaryContent::copyDict
            (
                input,
                allow,
                deny
            )
        );

        IOobject::writeDivider(Info);
        input.writeEntry("input", Info);
        copied.writeEntry("copied", Info);
    }

    for (label argi=1; argi < args.size(); ++argi)
    {
        const auto dictFile = args.get<fileName>(argi);
        IFstream is(dictFile);

        dictionary input(is);

        dictionary copied
        (
            dictionaryContent::copyDict
            (
                input,
                allow,
                deny
            )
        );

        IOobject::writeDivider(Info);
        input.writeEntry("input", Info);
        copied.writeEntry("copied", Info);
    }

    IOobject::writeEndDivider(Info);

    return 0;
}


// ************************************************************************* //
