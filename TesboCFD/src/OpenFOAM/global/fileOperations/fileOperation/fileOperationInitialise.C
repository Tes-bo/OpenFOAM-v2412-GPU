/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "fileOperationInitialise.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(fileOperationInitialise, 0);
    defineRunTimeSelectionTable(fileOperationInitialise, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::fileOperationInitialise::fileOperationInitialise
(
    int& argc,
    char**& argv
)
{
    // Check for -ioRanks: requires an argument
    int index = -1;
    for (int argi = 1; argi < argc; ++argi)
    {
        const char *optName = argv[argi];

        if (optName[0] == '-')
        {
            ++optName;  // Looks like an option, skip leading '-'
            bool emitErrorMessage = false;

            if (strcmp(optName, "ioRanks") == 0)
            {
                // Requires a parameter
                if (argi < argc-1)
                {
                    index = argi;
                    ++argi;
                    Foam::setEnv("FOAM_IORANKS", argv[argi], true);
                    break;
                }
                else
                {
                    emitErrorMessage = true;
                }
            }

            if (emitErrorMessage)
            {
                // Missing argument: emit message but not exit or
                // FatalError since Pstream etc are not yet initialised

                Info<< nl
                    << "Error: option '-" << optName
                    << "' requires an argument" << nl << nl;

                //NO: UPstream::exit(1);  // works for serial and parallel
            }
        }
    }

    if (index > 0)
    {
        for (int i=index+2; i<argc; i++)
        {
            argv[i-2] = argv[i];
        }
        argc -= 2;
    }
}


Foam::autoPtr<Foam::fileOperations::fileOperationInitialise>
Foam::fileOperations::fileOperationInitialise::New
(
    const word& type,
    int& argc,
    char**& argv
)
{
    DebugInFunction << "Constructing fileOperationInitialise" << endl;

    auto* ctorPtr = wordConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "fileOperationInitialise",
            type,
            *wordConstructorTablePtr_
        ) << abort(FatalError);
    }

    return autoPtr<fileOperationInitialise>(ctorPtr(argc, argv));
}


// ************************************************************************* //
