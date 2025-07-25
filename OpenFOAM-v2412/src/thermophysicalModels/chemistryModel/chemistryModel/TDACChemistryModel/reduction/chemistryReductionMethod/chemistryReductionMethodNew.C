/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "chemistryReductionMethod.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::chemistryReductionMethod<CompType, ThermoType>>
Foam::chemistryReductionMethod<CompType, ThermoType>::New
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
{
    const dictionary& reductionDict = dict.subDict("reduction");

    const word methodName(reductionDict.get<word>("method"));

    Info<< "Selecting chemistry reduction method " << methodName << endl;

    const word methodTypeName
    (
        methodName
      + '<' + CompType::typeName + ',' + ThermoType::typeName() + '>'
    );

    const auto& cnstrTable = *(dictionaryConstructorTablePtr_);

    auto* ctorPtr = cnstrTable.lookup(methodTypeName, nullptr);

    if (!ctorPtr)
    {
        const wordList names(cnstrTable.sortedToc());

        constexpr const int nCmpt = 7;

        /// DynamicList<word> thisCmpts(6);
        /// thisCmpts.append(CompType::typeName);
        /// thisCmpts.append
        /// (
        ///     basicThermo::splitThermoName(ThermoType::typeName(), 5)
        /// );
        ///
        /// DynamicList<word> validNames;

        DynamicList<wordList> validCmpts;
        validCmpts.append
        (
            // Header
            wordList
            ({
                typeName_(),
                "reactionThermo",
                "transport",
                "thermo",
                "equationOfState",
                "specie",
                "energy"
            })
        );

        for (const word& validName : names)
        {
            wordList cmpts(basicThermo::splitThermoName(validName, nCmpt));

            if (!cmpts.empty())
            {
                /// if (thisCmpts == SubList<word>(cmpts, 6, 1))
                /// {
                ///     validNames.append(cmpts[0]);
                /// }
                validCmpts.append(std::move(cmpts));
            }
        }

        FatalErrorInLookup
        (
            typeName_(),
            methodName,
            cnstrTable
        );

        if (validCmpts.size() > 1)
        {
            FatalError
                << "All " << validCmpts[0][0] << '/' << validCmpts[0][1]
                << "/thermoPhysics combinations:" << nl << nl;

            // Table of available packages (as constituent parts)
            printTable(validCmpts, FatalError) << nl;
        }

        FatalError
            << exit(FatalError);
    }

    return autoPtr<chemistryReductionMethod<CompType, ThermoType>>
    (
        ctorPtr(dict, chemistry)
    );
}


// ************************************************************************* //
