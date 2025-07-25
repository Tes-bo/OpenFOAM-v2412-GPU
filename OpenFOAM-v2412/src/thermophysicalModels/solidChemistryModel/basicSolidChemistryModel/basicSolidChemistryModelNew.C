/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "basicSolidChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSolidChemistryModel>
Foam::basicSolidChemistryModel::New(solidReactionThermo& thermo)
{
    const IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("chemistryType");

    Info<< "Selecting chemistry type " << chemistryTypeDict << endl;

    const IOdictionary thermoDict
    (
        IOobject
        (
            basicThermo::dictName,
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const dictionary& solidThermoTypeDict = thermoDict.subDict("thermoType");
    const word solidThermoTypeName
    (
        solidThermoTypeDict.get<word>("transport") + '<'
      + solidThermoTypeDict.get<word>("thermo") + '<'
      + solidThermoTypeDict.get<word>("equationOfState") + '<'
      + solidThermoTypeDict.get<word>("specie") + ">>,"
      + solidThermoTypeDict.get<word>("energy") + ">"
    );

    const dictionary& gasThermoTypeDict = thermoDict.subDict("gasThermoType");
    const word gasThermoTypeName
    (
        gasThermoTypeDict.get<word>("transport") + '<'
      + gasThermoTypeDict.get<word>("thermo") + '<'
      + gasThermoTypeDict.get<word>("equationOfState") + '<'
      + gasThermoTypeDict.get<word>("specie") + ">>,"
      + gasThermoTypeDict.get<word>("energy") + ">"
    );

    // Construct the name of the chemistry type from the components
    const word chemistryTypeName
    (
        chemistryTypeDict.get<word>("chemistrySolver") + '<'
      + chemistryTypeDict.get<word>("chemistryThermo") + '<'
      + typeName + ','
      + solidThermoTypeName + ',' + gasThermoTypeName + ">>"
    );

    Info<< "chemistryTypeName " << chemistryTypeName << endl;

    const auto& cnstrTable = *(thermoConstructorTablePtr_);

    auto* ctorPtr = cnstrTable.lookup(chemistryTypeName, nullptr);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            chemistryTypeDict,
            typeName,
            word::null, // Suppress long name? Just output dictionary (above)
            cnstrTable
        );

        // Table of available packages (as constituent parts)
        basicThermo::printThermoNames
        (
            FatalIOError,
            wordList
            ({
                "chemistrySolver",
                "chemistryThermo",
                "baseChemistry",
                "transport",
                "thermo",  // solid
                "equationOfState",
                "specie",
                "energy",
                "transport",
                "thermo",  // gas
                "equationOfState",
                "specie",
                "energy"
            }),
            cnstrTable.sortedToc()
        );

        FatalIOError
            << exit(FatalIOError);
    }

    return autoPtr<basicSolidChemistryModel>(ctorPtr(thermo));
}


// ************************************************************************* //
