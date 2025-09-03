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

#include "runTimeControl.H"
#include "dictionary.H"
#include "runTimeCondition.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(runTimeControl, 0);
    addToRunTimeSelectionTable(functionObject, runTimeControl, dictionary);
}
}
}

Foam::Enum
<
    Foam::functionObjects::runTimeControls::runTimeControl::satisfiedAction
>
Foam::functionObjects::runTimeControls::runTimeControl::satisfiedActionNames
{
    { satisfiedAction::ABORT, "abort"},
    { satisfiedAction::END, "end"},
    { satisfiedAction::SET_TRIGGER, "setTrigger"},
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeControl::runTimeControl
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    conditions_(),
    groupMap_(),
    nWriteStep_(0),
    writeStepI_(0),
    satisfiedAction_(satisfiedAction::END),
    triggerIndex_(labelMin),
    active_(getProperty("active", true)),
    canRestart_(getProperty("canRestart", false))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::runTimeControl::read
(
    const dictionary& dict
)
{
    if (functionObject::postProcess)
    {
        Info<< "Deactivated " << name()
            << " function object for post-processing"
            << endl;

        return false;
    }


    if (fvMeshFunctionObject::read(dict))
    {
        Info<< type() << " " << name() << ":" << nl;

        const dictionary& conditionsDict = dict.subDict("conditions");
        const wordList conditionNames(conditionsDict.toc());
        conditions_.setSize(conditionNames.size());

        label uniqueGroupi = 0;
        forAll(conditionNames, conditioni)
        {
            const word& conditionName = conditionNames[conditioni];
            const dictionary& dict = conditionsDict.subDict(conditionName);

            conditions_.set
            (
                conditioni,
                runTimeCondition::New(conditionName, obr_, dict, *this)
            );

            label groupi = conditions_[conditioni].groupID();

            if (groupMap_.insert(groupi, uniqueGroupi))
            {
                ++uniqueGroupi;
            }
        }

        dict.readIfPresent("nWriteStep", nWriteStep_);

        // Check that some conditions are set
        if (conditions_.empty())
        {
            Info<< "    No conditions present" << endl;
        }
        else
        {
            // Check that at least one condition is active
            bool check = false;
            for (const auto& condition : conditions_)
            {
                if (condition.active())
                {
                    check = true;
                    break;
                }
            }

            if (!check)
            {
                Info<< "    All conditions are inactive" << endl;
            }
        }

        Info<< endl;

        // Set the action to perform when all conditions are satisfied
        // - set to end for backwards compatibility with v1806
        satisfiedAction_ =
            satisfiedActionNames.getOrDefault
            (
                "satisfiedAction",
                dict,
                satisfiedAction::END
            );

        if (satisfiedAction_ == satisfiedAction::SET_TRIGGER)
        {
            triggerIndex_ = dict.get<label>("trigger");
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::runTimeControls::runTimeControl::execute()
{
    if (canRestart_)
    {
        active_ = true;
        canRestart_ = false;
    }

    if (!active_)
    {
        return true;
    }

    Info<< type() << " " << name() << " output:" << nl;

    // IDs of satisfied conditions
    DynamicList<label> IDs(conditions_.size());

    // Run stops only if all conditions within a group are satisfied
    List<bool> groupSatisfied(groupMap_.size(), true);
    List<bool> groupActive(groupMap_.size(), false);

    forAll(conditions_, conditioni)
    {
        runTimeCondition& condition = conditions_[conditioni];

        if (condition.active())
        {
            bool conditionSatisfied = condition.apply();

            const label groupi = condition.groupID();

            auto conditionIter = groupMap_.cfind(groupi);

            if (!conditionIter.good())
            {
                FatalErrorInFunction
                    << "group " << groupi << " not found in map"
                    << abort(FatalError);
            }

            if (conditionSatisfied)
            {
                IDs.append(conditioni);

                groupActive[conditionIter()] = true;

                if (groupi == -1)
                {
                    // Condition not part of a group - only requires this to be
                    // satisfied for completion flag to be set
                    groupSatisfied[conditionIter()] = true;
                    break;
                }
            }
            else
            {
                groupSatisfied[conditionIter()] = false;
            }
        }
    }

    bool done = false;
    forAll(groupSatisfied, groupi)
    {
        if (groupSatisfied[groupi] && groupActive[groupi])
        {
            done = true;
            break;
        }
    }

    if (done)
    {
        for (label conditioni : IDs)
        {
            Info<< "    " << conditions_[conditioni].type() << ": "
                << conditions_[conditioni].name()
                << " condition satisfied" << nl;
        }

        switch (satisfiedAction_)
        {
            case satisfiedAction::ABORT:
            case satisfiedAction::END:
            {
                // Set to write a data dump or finalise the calculation
                Time& time = const_cast<Time&>(time_);

                if (writeStepI_ < nWriteStep_ - 1)
                {
                    ++writeStepI_;
                    Info<< "    Writing fields - step " << writeStepI_ << nl;
                    time.writeNow();
                }
                else
                {
                    Info<< "    Stopping calculation" << nl
                        << "    Writing fields";

                    if (nWriteStep_ != 0)
                    {
                        Info<< " - final step";
                    }

                    Info<< nl << endl;
                    active_ = false;

                    // Write any registered objects and set the end-time
                    time.writeAndEnd();

                    // Trigger any function objects
                    time.run();

                    if (satisfiedAction_ == satisfiedAction::ABORT)
                    {
                        FatalErrorInFunction
                            << "Abort triggered"
                            << exit(FatalError);
                    }
                }
                break;
            }
            case satisfiedAction::SET_TRIGGER:
            {
                Info<< "    Setting trigger " << triggerIndex_ << nl;

                setTrigger(triggerIndex_);

                // Deactivate the model
                active_ = false;
                setProperty("active", active_);

                // Can be restarted
                canRestart_ = true;
                setProperty("canRestart", canRestart_);

                // Reset all conditions in case the control is recycled/trigger
                // index is set to a smaller value
                forAll(conditions_, conditioni)
                {
                    runTimeCondition& condition = conditions_[conditioni];
                    condition.reset();
                }

                break;
            }
        }
    }
    else
    {
        Info<< "    Conditions not met" << nl;
    }

    Info<< endl;

    return true;
}


bool Foam::functionObjects::runTimeControls::runTimeControl::write()
{
    for (auto& condition : conditions_)
    {
        condition.write();
    }

    return true;
}


// ************************************************************************* //
