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

#include "minMaxCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<>
void Foam::functionObjects::runTimeControls::minMaxCondition::
setValue<Foam::scalar>
(
    const word& valueType,
    const word& fieldName,
    scalar& value
) const
{
    state_.getObjectResult(functionObjectName_, fieldName, value);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(minMaxCondition, 0);
    addToRunTimeSelectionTable(runTimeCondition, minMaxCondition, dictionary);

}
}
}

const Foam::Enum
<
    Foam
  ::functionObjects
  ::runTimeControls
  ::minMaxCondition
  ::modeType
>
Foam::functionObjects::runTimeControls::minMaxCondition::modeTypeNames_
({
    { modeType::mdMin, "minimum" },
    { modeType::mdMax, "maximum" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::minMaxCondition::minMaxCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    functionObjectName_(dict.get<word>("functionObject")),
    mode_(modeTypeNames_.get("mode", dict)),
    fieldNames_(dict.get<wordList>("fields")),
    value_(dict.get<scalar>("value"))
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::minMaxCondition::apply()
{
    bool satisfied = true;

    if (!active_)
    {
        return satisfied;
    }

    for (const word& fieldName :fieldNames_)
    {
        const word valueType =
            state_.objectResultType(functionObjectName_, fieldName);

        if (valueType.empty())
        {
            WarningInFunction
                << "Unable to find entry " << fieldName
                << " for function object " << functionObjectName_
                << ".  Condition will not be applied."
                << endl;

            continue;
        }

        scalar v = 0;
        setValue<scalar>(valueType, fieldName, v);
        setValue<vector>(valueType, fieldName, v);
        setValue<sphericalTensor>(valueType, fieldName, v);
        setValue<symmTensor>(valueType, fieldName, v);
        setValue<tensor>(valueType, fieldName, v);

        Switch ok = false;
        switch (mode_)
        {
            case mdMin:
            {
                if (v < value_)
                {
                    ok = true;
                }
                break;
            }
            case mdMax:
            {
                if (v > value_)
                {
                    ok = true;
                }
                break;
            }
        }

        Log << "    " << type() << ": " << modeTypeNames_[mode_] << " "
            << fieldName << ": value = " << v
            << ", threshold value = " << value_
            << ", satisfied = " << ok << endl;

        satisfied = satisfied && ok;
    }

    return satisfied;
}


void Foam::functionObjects::runTimeControls::minMaxCondition::write()
{
    // do nothing
}


void Foam::functionObjects::runTimeControls::minMaxCondition::reset()
{
    // do nothing
}


// ************************************************************************* //
