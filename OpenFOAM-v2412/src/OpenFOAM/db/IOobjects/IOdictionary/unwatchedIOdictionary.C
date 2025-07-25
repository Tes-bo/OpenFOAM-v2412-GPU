/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "unwatchedIOdictionary.H"
#include "objectRegistry.H"
#include "Pstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unwatchedIOdictionary::unwatchedIOdictionary
(
    const IOobject& io,
    const dictionary* fallback
)
:
    unwatchedIOdictionary(io, typeName, fallback)
{}


Foam::unwatchedIOdictionary::unwatchedIOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    unwatchedIOdictionary(io, typeName, &dict)
{}


Foam::unwatchedIOdictionary::unwatchedIOdictionary
(
    const IOobject& io,
    const word& wantedType,
    const dictionary* fallback
)
:
    baseIOdictionary(io, fallback)
{
    if (!readHeaderOk(IOstreamOption::ASCII, wantedType) && fallback)
    {
        dictionary::operator=(*fallback);
    }

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


Foam::unwatchedIOdictionary::unwatchedIOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    baseIOdictionary(io, is)
{
    // Default construct dictionary and read in afterwards
    // so that if there is some fancy massaging due to a
    // functionEntry in
    // the dictionary at least the type information is already complete.
    is  >> *this;

    // For if MUST_READ_IF_MODIFIED
    addWatch();
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

Foam::label Foam::unwatchedIOdictionary::addWatch(const fileName& f)
{
    label index = -1;

    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        index = files_.find(f);

        if (index == -1)
        {
            index = files_.size();
            files_.append(f);
        }
    }
    return index;
}


void Foam::unwatchedIOdictionary::addWatch()
{
    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        fileName f = filePath();
        if (f.empty())
        {
            // We don't have this file but would like to re-read it.
            // Possibly if master-only reading mode.
            f = objectPath();
        }

        if (files_.found(f))
        {
            FatalErrorInFunction
                << "Object " << objectPath() << " of type " << type()
                << " already watched" << nl
                << abort(FatalError);
        }

        // If master-only reading only the master will have all dependencies
        // so broadcast these to other ranks
        const bool masterOnly
        (
            global()
         && (
                IOobject::fileModificationChecking == IOobject::timeStampMaster
             || IOobject::fileModificationChecking == IOobject::inotifyMaster
            )
        );

        if (masterOnly && Pstream::parRun())
        {
            Pstream::broadcast(files_);
        }

        addWatch(f);
    }
}


// ************************************************************************* //
