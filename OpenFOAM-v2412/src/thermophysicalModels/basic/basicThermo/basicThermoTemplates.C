/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "basicThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class ThermoConstructTable>
typename ThermoConstructTable::mapped_type
Foam::basicThermo::getThermoOrDie
(
    const dictionary& thermoTypeDict,
    ThermoConstructTable& thermoTable,
    const word& thermoTypeName,
    const wordList& cmptNames
)
{
    // Lookup the thermo package

    auto ctorIter = thermoTable.cfind(thermoTypeName);

    // Print error message if package not found in the table
    if (!ctorIter.good())
    {
        FatalIOErrorInLookup
        (
            thermoTypeDict,
            Thermo::typeName,
            word::null, // Suppress long name? Just output dictionary (above)
            thermoTable
        );

        basicThermo::printThermoNames
        (
            FatalIOError,
            cmptNames,
            thermoTable.sortedToc()
        ) << exit(FatalIOError);

        // return nullptr;
    }

    return ctorIter.val();
}


template<class Thermo, class ThermoConstructTable>
typename ThermoConstructTable::mapped_type
Foam::basicThermo::getThermoOrDie
(
    const dictionary& thermoDict,
    ThermoConstructTable& thermoTable
)
{
    const dictionary* dictptr = thermoDict.findDict("thermoType");

    if (dictptr)
    {
        const auto& thermoTypeDict = *dictptr;

        const wordList* cmptHeaderPtr = &(wordList::null());

        // Thermo package name, constructed from components
        const word thermoTypeName
        (
            basicThermo::makeThermoName(thermoTypeDict, cmptHeaderPtr)
        );

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        return getThermoOrDie<Thermo, ThermoConstructTable>
        (
            thermoTypeDict,
            thermoTable,
            thermoTypeName,
            *cmptHeaderPtr
        );
    }
    else
    {
        const word thermoTypeName(thermoDict.get<word>("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        auto ctorIter = thermoTable.cfind(thermoTypeName);

        if (!ctorIter.good())
        {
            FatalIOErrorInLookup
            (
                thermoDict,
                Thermo::typeName,
                thermoTypeName,
                thermoTable
            ) << exit(FatalIOError);
        }

        return ctorIter.val();
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto* ctorPtr = getThermoOrDie<Thermo>
    (
        thermoDict,
        *(Thermo::fvMeshConstructorTablePtr_)
    );

    return autoPtr<Thermo>(ctorPtr(mesh, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    auto* ctorPtr = getThermoOrDie<Thermo>
    (
        dict,
        *(Thermo::dictionaryConstructorTablePtr_)
    );

    return autoPtr<Thermo>(ctorPtr(mesh, dict, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto* ctorPtr = getThermoOrDie<Thermo>
    (
        thermoDict,
        *(Thermo::fvMeshDictPhaseConstructorTablePtr_)
    );

    return autoPtr<Thermo>(ctorPtr(mesh, phaseName, dictName));
}


// ************************************************************************* //
