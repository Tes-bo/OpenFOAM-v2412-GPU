/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "valueAverageBase.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::functionObjects::valueAverageBase::windowType
>
Foam::functionObjects::valueAverageBase::windowTypeNames
({
    { windowType::NONE, "none" },
    { windowType::APPROXIMATE, "approximate" },
    { windowType::EXACT, "exact" }
});


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

void Foam::functionObjects::valueAverageBase::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Value averages");
    writeCommented(os, "Time");
    forAll(fieldNames_, fieldi)
    {
        writeTabbed(os, fieldNames_[fieldi]);
    }
    os  << endl;
}


void Foam::functionObjects::valueAverageBase::readState(dictionary& dict)
{
    if (resetOnRestart_)
    {
        resetState(dict);
        return;
    }

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName = fieldNames_[fieldi];

        if (dict.found(fieldName))
        {
            const dictionary& valueDict = dict.subDict(fieldName);
            valueDict.readEntry("totalTime", totalTime_[fieldi]);
        }
        else
        {
            dict.set(fieldName, dictionary());
            totalTime_[fieldi] = 0;
        }
    }
}


void Foam::functionObjects::valueAverageBase::writeState(dictionary& dict)
{
    forAll(fieldNames_, fieldi)
    {
        const word& fieldName = fieldNames_[fieldi];

        if (dict.found(fieldName))
        {
            dictionary& valueDict = dict.subDict(fieldName);
            valueDict.add("totalTime", totalTime_[fieldi], true);
        }
        else
        {
            dictionary valueDict;
            valueDict.add("totalTime", totalTime_[fieldi], true);
            dict.add(fieldName, valueDict);
        }
    }
}


void Foam::functionObjects::valueAverageBase::resetState(dictionary& dict)
{
    forAll(fieldNames_, fieldi)
    {
        dict.set(fieldNames_[fieldi], dictionary());
        totalTime_[fieldi] = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::valueAverageBase::valueAverageBase
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state,
    const bool writeToFile
)
:
    writeFile(obr, name, state.type(), dict, writeToFile),
    log(true),
    resetOnRestart_(false),
    windowType_(windowType::NONE),
    state_(state),
    functionObjectName_("unknown-functionObject"),
    fieldNames_(),
    tolerance_(dict.getOrDefault<scalar>("tolerance", -1)),
    window_(-1),
    totalTime_()
{
    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::valueAverageBase::read(const dictionary& dict)
{
    if (writeFile::read(dict))
    {
        // Make certain that the values are consistent with the defaults:
        resetOnRestart_ = false;

        dict.readEntry("functionObject", functionObjectName_);
        dict.readEntry("fields", fieldNames_);
        if (dict.readIfPresent("window", window_))
        {
            window_ = state_.time().userTimeToTime(window_);

            if (window_ > 0)
            {
                windowType_ = windowTypeNames.get("windowType", dict);
            }
        }

        totalTime_.resize(fieldNames_.size(), Zero);

        dict.readIfPresent("resetOnRestart", resetOnRestart_);

        dict.readIfPresent("log", log);

        return true;
    }

    return false;
}


bool Foam::functionObjects::valueAverageBase::calculate(dictionary& dict)
{
    scalar dt = state_.time().deltaTValue();

    Log << indent << state_.type() << ": " << prefix_.c_str()
        << " averages:" << nl;

    file() << state_.time().timeName();

    DynamicList<word> unprocessedFields(fieldNames_.size());

    bool converged = true;

    forAll(fieldNames_, fieldi)
    {
        totalTime_[fieldi] += dt;

        const bool processed =
        (
             calc<label, scalar>(fieldi, converged, dict)
          || calc<scalar>(fieldi, converged, dict)
          || calc<vector>(fieldi, converged, dict)
          || calc<sphericalTensor>(fieldi, converged, dict)
          || calc<symmTensor>(fieldi, converged, dict)
          || calc<tensor>(fieldi, converged, dict)
        );

        if (!processed)
        {
            unprocessedFields.append(fieldNames_[fieldi]);

            file() << tab << "n/a";
        }
    }

    file() << endl;

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        for (const word& fieldName : unprocessedFields)
        {
            Log << indent << "        " << fieldName << nl;
        }

        if (unprocessedFields.size() == fieldNames_.size())
        {
            converged = false;
        }
    }

    Log << endl;

    return converged;
}


// ************************************************************************* //
