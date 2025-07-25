/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "zeroATCcells.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(zeroATCcells, 0);
defineRunTimeSelectionTable(zeroATCcells, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zeroATCcells::zeroATCcells
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    zeroATCPatches_(),
    zeroATCZones_(),
    zeroATCcells_()
{
    dict.readIfPresent("zeroATCPatchTypes", zeroATCPatches_);

    wordList zeroATCZoneNames;

    if (dict.readIfPresent("zeroATCZones", zeroATCZoneNames))
    {
        zeroATCZones_.resize(zeroATCZoneNames.size(), -1);

        forAll(zeroATCZoneNames, zI)
        {
            const word& zoneName = zeroATCZoneNames[zI];

            label zoneID = mesh.cellZones().findZoneID(zoneName);
            if (zoneID == -1)
            {
                WarningInFunction
                    << "cannot find cellZone "
                    << zoneName
                    << " for smoothing ATC"
                    << endl;
            }
            zeroATCZones_[zI] = zoneID;
        }
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<zeroATCcells> zeroATCcells::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType
    (
        dict.getOrDefault<word>("maskType", "faceCells")
    );

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "zeroATCcells",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<zeroATCcells> (ctorPtr(mesh,dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
