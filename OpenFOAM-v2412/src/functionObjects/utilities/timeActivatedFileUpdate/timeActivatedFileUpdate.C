/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "timeActivatedFileUpdate.H"
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeActivatedFileUpdate, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        timeActivatedFileUpdate,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::timeActivatedFileUpdate::updateFile()
{
    modified_ = false;

    label i = lastIndex_;
    while
    (
        i < timeVsFile_.size()-1
     && timeVsFile_[i+1].first() < time_.value()+0.5*time_.deltaTValue()
    )
    {
        i++;
    }

    if (i > lastIndex_)
    {
        const fileName& srcFile = timeVsFile_[i].second();

        // Report case-relative path for information
        Log << nl << type() << ": copying file" << nl
            << "from: " << time_.relativePath(srcFile, true) << nl
            << "to  : " << time_.relativePath(fileToUpdate_, true) << nl
            << endl;

        if
        (
            UPstream::master()
         || (
                fileHandler().distributed()
             && UPstream::master(fileHandler().comm())
            )
        )
        {
            // Copy on master only for non-distributed
            fileName tmpFile(fileToUpdate_ + Foam::name(pid()));
            Foam::cp(srcFile, tmpFile);
            Foam::mv(tmpFile, fileToUpdate_);
        }
        lastIndex_ = i;
        modified_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeActivatedFileUpdate::timeActivatedFileUpdate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    fileToUpdate_(),
    timeVsFile_(),
    lastIndex_(-1),
    modified_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeActivatedFileUpdate::read
(
    const dictionary& dict
)
{
    timeFunctionObject::read(dict);

    dict.readEntry("fileToUpdate", fileToUpdate_);
    dict.readEntry("timeVsFile", timeVsFile_);

    lastIndex_ = -1;
    fileToUpdate_.expand();

    if (fileToUpdate_.empty() || timeVsFile_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Bad entries for fileToUpdate and/or timeVsFile" << endl
            << exit(FatalIOError);
    }

    Info<< type() << " " << name() << " output:" << nl
        << "    time vs file list:" << nl;

    for (auto& tuple : timeVsFile_)
    {
        fileName& srcFile = tuple.second();
        srcFile.expand();

        // Report case-relative path for information
        Info<< "    " << tuple.first() << tab
            << time_.relativePath(srcFile, true) << nl;

        // No need for distributed test since dictionaries are read on the
        // master only or otherwise they need to be copied everywhere
        if (UPstream::master())  // || time_.distributed())
        {
            if (!Foam::isFile(srcFile))
            {
                // Report full path on error
                FatalIOErrorInFunction(dict)
                    << "File not found: " << srcFile << endl
                    << exit(FatalIOError);
            }
        }
    }

    // Copy starting files
    updateFile();

    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::execute()
{
    updateFile();

    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::write()
{
    return true;
}


bool Foam::functionObjects::timeActivatedFileUpdate::filesModified() const
{
    return modified_;
}


// ************************************************************************* //
