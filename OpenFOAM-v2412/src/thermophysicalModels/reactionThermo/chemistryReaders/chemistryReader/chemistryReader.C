/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "chemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::chemistryReader<ThermoType>>
Foam::chemistryReader<ThermoType>::New
(
    const dictionary& thermoDict,
    speciesTable& species
)
{
    // Use specified reader or default to CHEMKIN for backward compatibility
    const word readerName
    (
        thermoDict.getOrDefault<word>("chemistryReader", "chemkinReader")
    );

    Info<< "Selecting chemistryReader " << readerName << endl;

    auto* ctorPtr = dictionaryConstructorTable(readerName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            thermoDict,
            "chemistryReader",
            readerName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<chemistryReader<ThermoType>>
    (
        ctorPtr(thermoDict, species)
    );
}


// ************************************************************************* //
