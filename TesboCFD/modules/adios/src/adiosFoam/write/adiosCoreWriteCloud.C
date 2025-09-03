/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "adiosCoreWrite.H"

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

Foam::label Foam::adiosFoam::adiosCoreWrite::writeCloudRegistry
(
    const cloud* cldPtr,
    objectRegistry& obr
)
{
    obr.clear();

    if (!cldPtr)
    {
        return 0;
    }

    cldPtr->writeObjects(obr);

    const auto* pointsPtr = cloud::findIOPosition(obr);

    if (!pointsPtr)
    {
        // This should be impossible
        WarningInFunction
            << "Cloud " << cldPtr->name() << " did not write 'position'"
            << nl;
        return 0;
    }

    const label nParcels = returnReduce(cldPtr->nParcels(), sumOp<label>());

    const adiosFoam::cloudInfo cInfo(*cldPtr);
    const fileName varName(cInfo.fullName());


    // Cloud attributes:
    // - class
    // - number of parcels (all processes) as an attribute
    //
    if (UPstream::master())
    {
        putAttribute(varName/"class",         cInfo.type());
        putIntAttribute(varName/"nParcels",   nParcels);
    }

    DynamicList<word> written(obr.size());

    written.append(writeCloudFields<label>(varName, obr));
    written.append(writeCloudFields<scalar>(varName, obr));
    written.append(writeCloudFields<vector>(varName, obr));
    written.append(writeCloudFields<symmTensor>(varName, obr));
    written.append(writeCloudFields<sphericalTensor>(varName, obr));
    written.append(writeCloudFields<tensor>(varName, obr));

    obr.erase(written);

    if (obr.size())
    {
        Info<< "Did not write all cloud fields to ADIOS: "
            << flatOutput(obr.sortedNames()) << nl;
    }

    return written.size();
}


// ************************************************************************* //
