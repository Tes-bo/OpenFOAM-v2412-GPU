/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "abort.H"
#include "dictionary.H"
#include "error.H"
#include "Time.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(abort, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        abort,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read file contents and return a stop control as follows:
//
// - action=writeNow, action=nextWrite action=noWriteNow :
//   The signalled action. Report as corresponding <action>.
//
// Anything else (empty file, no action=, etc) is reported as <unknown>.
//
static enum Time::stopAtControls getStopAction(const std::string& filename)
{
    // Slurp entire input file (must exist) as a single string
    std::string fileContent;

    std::ifstream is(filename);
    std::getline(is, fileContent, '\0');

    const auto equals = fileContent.find('=');

    if (equals != std::string::npos)
    {
        const word actionName(word::validate(fileContent.substr(equals+1)));

        return
            Time::stopAtControlNames.lookup
            (
                actionName,
                Time::stopAtControls::saUnknown
            );
    }

    return Time::stopAtControls::saUnknown;
}


// Long description for the action name
static std::string longDescription(const Time::stopAtControls ctrl)
{
    switch (ctrl)
    {
        case Time::saEndTime :
        {
            return "continue simulation to the endTime";
            break;
        }

        case Time::saNoWriteNow :
        {
            return "stop without writing data";
            break;
        }

        case Time::saWriteNow :
        {
            return "stop and write data";
            break;
        }

        case Time::saNextWrite :
        {
            return "stop after next data write";
            break;
        }

        default:
        {
            // Invalid choices already filtered out by Enum
            return "unknown action";
            break;
        }
    }
}

}  // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::abort::abort
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    file_(),
    defaultAction_(Time::stopAtControls::saUnknown),
    triggered_(false)
{
    read(dict);

    // Cleanup old files from previous runs
    if (Pstream::master())
    {
        Foam::rm(file_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::abort::read(const dictionary& dict)
{
    timeFunctionObject::read(dict);

    file_.clear();

    if (dict.readIfPresent("file", file_))
    {
        file_.expand();
        if (!file_.empty() && !file_.isAbsolute())
        {
            file_ = time_.globalPath()/file_;
            file_.clean();  // Remove unneeded ".."
        }
    }

    // Ensure we always have a reasonable default file
    if (file_.empty())
    {
        file_ = time_.globalPath()/name();
        file_.clean();  // Remove unneeded ".."
    }

    triggered_ = false;

    defaultAction_ = Time::stopAtControlNames.getOrDefault
    (
        "action",
        dict,
        Time::stopAtControls::saNextWrite
    );

    Info<< type() << " activated ("
        << longDescription(defaultAction_).c_str() <<")" << nl
        << "    File: " << file_ << endl;

    return true;
}


bool Foam::functionObjects::abort::execute()
{
    // If it has been triggered (eg, nextWrite) don't need to check it again
    if (!triggered_)
    {
        auto action = Time::stopAtControls::saUnknown;

        if (Pstream::master() && Foam::isFile(file_))
        {
            action = getStopAction(file_);

            if (Time::stopAtControls::saUnknown == action)
            {
                // An unknown action means an empty file or bad content.
                // Treat as a request for the default action.

                action = defaultAction_;
            }
        }

        // Send to sub-ranks. Also acts as an MPI barrier
        int intAction(action);
        Pstream::broadcast(intAction);

        action = Time::stopAtControls(intAction);

        // Call stopAt() on all processes
        triggered_ = time_.stopAt(action);

        if (triggered_)
        {
            Info<< "USER REQUESTED ABORT (timeIndex="
                << time_.timeIndex() << "): "
                << longDescription(action).c_str() << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::abort::write()
{
    return true;
}


bool Foam::functionObjects::abort::end()
{
    // Cleanup trigger file
    if (Pstream::master())
    {
        Foam::rm(file_);
    }

    return true;
}


// ************************************************************************* //
