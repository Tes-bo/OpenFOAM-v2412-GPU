/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "equationInitialResidualCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

#define SetResidual(Type)                                                     \
    setResidual<Type>(mesh, solverDict, fieldName, component, canSet, residual);

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(equationInitialResidualCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        equationInitialResidualCondition,
        dictionary
    );

}
}
}

const Foam::Enum
<
    Foam
  ::functionObjects
  ::runTimeControls
  ::equationInitialResidualCondition
  ::operatingMode
>
Foam::functionObjects::runTimeControls::equationInitialResidualCondition::
operatingModeNames
({
    { operatingMode::omMin, "minimum" },
    { operatingMode::omMax, "maximum" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::
equationInitialResidualCondition::equationInitialResidualCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    fieldSelection_(obr, true),
    value_(dict.get<scalar>("value")),
    timeStart_(dict.getOrDefault("timeStart", -GREAT)),
    mode_(operatingModeNames.get("mode", dict))
{
    fieldSelection_.read(dict);

    if (fieldSelection_.size())
    {
        timeStart_ = obr.time().userTimeToTime(timeStart_);
    }
    else
    {
        WarningInFunction
            << "No fields supplied: deactivating" << endl;

        active_ = false;
    }
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::
equationInitialResidualCondition::apply()
{
    fieldSelection_.updateSelection();

    bool satisfied = false;

    if (!active_)
    {
        return true;
    }

    if ((obr_.time().timeIndex() < 3) || (obr_.time().value() < timeStart_))
    {
        // Do not start checking until reached start time
        return false;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const dictionary& solverDict = mesh.data().solverPerformanceDict();

    const auto& selection = fieldSelection_.selection();
    List<scalar> result(selection.size(), -VGREAT);

    forAll(selection, fieldi)
    {
        const auto& fieldInfo = selection[fieldi];
        const word& fieldName = fieldInfo.name();

        if (solverDict.found(fieldName))
        {
            label component = fieldInfo.component();
            scalar residual = VGREAT;
            bool canSet = true;
            SetResidual(scalar);
            SetResidual(vector);
            SetResidual(symmTensor);
            SetResidual(sphericalTensor);
            SetResidual(tensor);

            result[fieldi] = residual;

            switch (mode_)
            {
                case omMin:
                {
                    if (residual < value_)
                    {
                        satisfied = true;
                    }
                    break;
                }
                case omMax:
                {
                    if (residual > value_)
                    {
                        satisfied = true;
                    }
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unhandled enumeration "
                        << operatingModeNames[mode_]
                        << abort(FatalError);
                }
            }
        }
    }

    bool valid = false;
    forAll(result, i)
    {
        if (result[i] < 0)
        {
            WarningInFunction
                << "Initial residual data not found for field "
                << selection[i].name()
                << ".  Solver dictionary contains " << solverDict.sortedToc()
                << endl;
        }
        else
        {
            valid = true;
        }
    }

    if (!valid)
    {
        WarningInFunction
            << "Initial residual data not found for any fields: "
            << "deactivating" << endl;

        active_ = false;
    }

    if (satisfied && valid)
    {
        Log << type() << ": " << name_
            << ": satisfied using threshold value: " << value_ << nl;

        forAll(result, resulti)
        {
            if (result[resulti] > 0)
            {
                Log << "    field: " << selection[resulti].name()
                    << ", residual: " << result[resulti] << nl;
            }
        }
        Log << endl;
    }

    return satisfied;
}


void Foam::functionObjects::runTimeControls::
equationInitialResidualCondition::write()
{
    // do nothing
}


void Foam::functionObjects::runTimeControls::
equationInitialResidualCondition::reset()
{
    // do nothing
}


// ************************************************************************* //
