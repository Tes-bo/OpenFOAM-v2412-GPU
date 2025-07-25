/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "systemCall.H"
#include "Time.H"
#include "dynamicCode.H"
#include "foamVersion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(systemCall, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        systemCall,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::systemCall::dispatch(const stringList& calls)
{
    if (calls.empty())
    {
        return 0;
    }

    label nCalls = 0;

    if (!masterOnly_ || Pstream::master())
    {
        for (const string& call : calls)
        {
            Foam::system(call); // Handles empty command as a successful no-op.
            ++nCalls;
        }
    }

    // MPI barrier
    if (masterOnly_)
    {
        Pstream::broadcast(nCalls);
    }

    return nCalls;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::systemCall::systemCall
(
    const word& name,
    const Time&,
    const dictionary& dict
)
:
    functionObject(name),
    executeCalls_(),
    writeCalls_(),
    endCalls_(),
    masterOnly_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::systemCall::read(const dictionary& dict)
{
    functionObject::read(dict);

    executeCalls_.clear();
    writeCalls_.clear();
    endCalls_.clear();

    dict.readIfPresent("executeCalls", executeCalls_);
    dict.readIfPresent("writeCalls", writeCalls_);
    dict.readIfPresent("endCalls", endCalls_);
    masterOnly_ = dict.getOrDefault("master", false);

    if (executeCalls_.empty() && endCalls_.empty() && writeCalls_.empty())
    {
        WarningInFunction
            << "No executeCalls, endCalls or writeCalls defined."
            << endl;
    }
    else if (isAdministrator())
    {
        FatalErrorInFunction
            << "System calls should not be executed by someone"
            << " with administrator rights for security reasons." << nl
            << nl << endl
            << exit(FatalError);
    }
    else if (!dynamicCode::allowSystemOperations)
    {
        FatalErrorInFunction
            << "Executing user-supplied system calls may have been disabled"
            << " by default" << nl
            << "for security reasons." << nl
            << "If you trust the code, you may enable this by adding"
            << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "to the InfoSwitches setting in the system controlDict." << nl
            << "The system controlDict is any of" << nl << nl
            << "    ~/.OpenFOAM/" << foamVersion::api << "/controlDict" << nl
            << "    ~/.OpenFOAM/controlDict" << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl << endl
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::systemCall::execute()
{
    dispatch(executeCalls_);
    return true;
}


bool Foam::functionObjects::systemCall::write()
{
    dispatch(writeCalls_);
    return true;
}


bool Foam::functionObjects::systemCall::end()
{
    dispatch(endCalls_);
    return true;
}


// ************************************************************************* //
