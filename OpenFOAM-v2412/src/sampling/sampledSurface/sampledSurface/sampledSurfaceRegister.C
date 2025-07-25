/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "sampledSurface.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polySurface* Foam::sampledSurface::getRegistrySurface
(
    const objectRegistry& obr,
    word lookupName
) const
{
    return
    (
        lookupName.empty()
      ? obr.getObjectPtr<polySurface>(this->name())
      : obr.getObjectPtr<polySurface>(lookupName)
    );
}


Foam::polySurface* Foam::sampledSurface::storeRegistrySurface
(
    objectRegistry& obr,
    word lookupName
) const
{
    if (lookupName.empty())
    {
        lookupName = this->name();
    }

    auto* surfptr = obr.getObjectPtr<polySurface>(lookupName);

    if (!surfptr)
    {
        surfptr = new polySurface(lookupName, obr);
        regIOobject::store(surfptr);
    }

    // Copy in geometry (removes existing fields if sizes have changed)
    surfptr->copySurface(*this);

    return surfptr;
}


bool Foam::sampledSurface::removeRegistrySurface
(
    objectRegistry& obr,
    word lookupName
) const
{
    return
    (
        lookupName.empty()
      ? polySurface::Delete(this->name(), obr)
      : polySurface::Delete(lookupName, obr)
    );
}


// ************************************************************************* //
