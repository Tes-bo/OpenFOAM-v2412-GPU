/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "pointToPoint.H"
#include "polyMesh.H"
#include "pointSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToPoint, 0);
    addToRunTimeSelectionTable(topoSetSource, pointToPoint, word);
    addToRunTimeSelectionTable(topoSetSource, pointToPoint, istream);
    addToRunTimeSelectionTable(topoSetPointSource, pointToPoint, word);
    addToRunTimeSelectionTable(topoSetPointSource, pointToPoint, istream);
}


Foam::topoSetSource::addToUsageTable Foam::pointToPoint::usage_
(
    pointToPoint::typeName,
    "\n    Usage: pointToPoint <pointSet>\n\n"
    "    Select all points in the pointSet\n\n"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointToPoint::pointToPoint
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSetPointSource(mesh),
    names_(Foam::one{}, setName),
    isZone_(false)
{}


Foam::pointToPoint::pointToPoint
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetPointSource(mesh, dict),
    names_(),
    isZone_(topoSetSource::readNames(dict, names_))
{}


Foam::pointToPoint::pointToPoint
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetPointSource(mesh),
    names_(Foam::one{}, word(checkIs(is))),
    isZone_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToPoint::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all elements of point "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            if (isZone_)
            {
                set.addSet(mesh_.pointZones()[setName]);
            }
            else
            {
                pointSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);

                set.addSet(loadedSet);
            }
        }
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all elements of point "
                << (isZone_ ? "zones: " : "sets: ")
                << flatOutput(names_) << nl;
        }

        for (const word& setName : names_)
        {
            if (isZone_)
            {
                set.subtractSet(mesh_.pointZones()[setName]);
            }
            else
            {
                pointSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);

                set.subtractSet(loadedSet);
            }
        }
    }
}


// ************************************************************************* //
