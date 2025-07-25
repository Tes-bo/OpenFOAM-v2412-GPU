/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
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

#include "regionProperties.H"
#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionProperties::regionProperties(const Time& runTime)
:
    regionProperties(runTime, IOobject::MUST_READ_IF_MODIFIED)
{}


Foam::regionProperties::regionProperties
(
    const Time& runTime,
    IOobjectOption::readOption rOpt
)
{
    HashTable<wordList>& props = *this;

    IOdictionary iodict
    (
        IOobject
        (
            "regionProperties",
            runTime.time().constant(),
            runTime.db(),
            rOpt,
            IOobjectOption::NO_WRITE
        )
    );

    if (IOobjectOption::isReadRequired(rOpt) || iodict.size())
    {
        iodict.readEntry("regions", props);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::regionProperties::count() const
{
    label n = 0;

    const HashTable<wordList>& props = *this;

    forAllConstIters(props, iter)
    {
        n += iter.val().size();
    }

    return n;
}


Foam::wordList Foam::regionProperties::names() const
{
    wordList list(this->count());

    label n = 0;

    const HashTable<wordList>& props = *this;

    for (const word& grp : props.sortedToc())
    {
        for (const word& name : props[grp])
        {
            list[n] = name;
            ++n;
        }
    }

    return list;
}


Foam::wordList Foam::regionProperties::sortedNames() const
{
    wordList list(this->count());

    label n = 0;

    const HashTable<wordList>& props = *this;

    forAllConstIters(props, iter)
    {
        for (const word& name : iter.val())
        {
            list[n] = name;
            ++n;
        }
    }

    Foam::sort(list);

    return list;
}


// ************************************************************************* //
