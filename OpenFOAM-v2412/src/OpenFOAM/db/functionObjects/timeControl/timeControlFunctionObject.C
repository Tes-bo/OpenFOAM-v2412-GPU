/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "timeControlFunctionObject.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeControl, 0);
}
}

const Foam::Enum<Foam::functionObjects::timeControl::controlMode>
Foam::functionObjects::timeControl::controlModeNames_
({
    { controlMode::TIME, "time" },
    { controlMode::TRIGGER, "trigger" },
    { controlMode::TIME_OR_TRIGGER, "timeOrTrigger" },
    { controlMode::TIME_AND_TRIGGER, "timeAndTrigger" }
});


// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::functionObjects::timeControl::readControls()
{
    if (dict_.readIfPresent("timeStart", timeStart_))
    {
        timeStart_ = time_.userTimeToTime(timeStart_);
    }
    if (dict_.readIfPresent("timeEnd", timeEnd_))
    {
        timeEnd_ = time_.userTimeToTime(timeEnd_);
    }

    controlMode_ =
        controlModeNames_.getOrDefault("controlMode", dict_, controlMode::TIME);
    dict_.readIfPresent("triggerStart", triggerStart_);
    dict_.readIfPresent("triggerEnd", triggerEnd_);

    deltaTCoeff_ = GREAT;
    if (dict_.readIfPresent("deltaTCoeff", deltaTCoeff_))
    {
        nStepsToStartTimeChange_ = labelMax;
    }
    else
    {
        nStepsToStartTimeChange_ = 3;
        dict_.readIfPresent
        (
            "nStepsToStartTimeChange",
            nStepsToStartTimeChange_
        );
    }
}


bool Foam::functionObjects::timeControl::active() const
{
    label triggeri = time_.functionObjects().triggerIndex();

    bool inTime =
        time_.value() >= (timeStart_ - 0.5*time_.deltaTValue())
     && time_.value() <= (timeEnd_ + 0.5*time_.deltaTValue());

    bool inTrigger = triggeri >= triggerStart_ && triggeri <= triggerEnd_;

    DebugInFunction
        << name() << " mode:" << controlModeNames_[controlMode_] << nl
        << " - time:" << time_.value()
        << " timeStart:" << timeStart_
        << " timeEnd:" << timeEnd_
        << " inTime:" << inTime << nl
        << " - triggeri:" << triggeri
        << " triggerStart:" << triggerStart_
        << " triggerEnd:" << triggerEnd_
        << " inTrigger:" << inTrigger
        << endl;

    switch (controlMode_)
    {
        case controlMode::TIME:
        {
            return inTime;
        }
        case controlMode::TRIGGER:
        {
            return inTrigger;
        }
        case controlMode::TIME_OR_TRIGGER:
        {
            return inTime || inTrigger;
        }
        case controlMode::TIME_AND_TRIGGER:
        {
            return inTime && inTrigger;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration: " << controlModeNames_[controlMode_]
                << abort(FatalError);
        }
    }

    return false;
}


Foam::scalar Foam::functionObjects::timeControl::calcExpansion
(
    const scalar startRatio,
    const scalar y,
    const label n
)
{
    scalar ratio = startRatio;
    // This function is used to calculate the 'expansion ratio' or
    // 'time step growth factor'. Inputs:
    //  - ratio: wanted ratio
    //  - y: overall time required (divided by deltaT) to get from 1 to ratio
    //  - n: number of steps in which to get to wanted ratio

    // Let a0 = first term in geometric progression (GP), an = last term,
    // n = number of terms,  ratio = ratio of GP, Sn = sum of series for
    // first n terms.
    // We have an=a0*ratio^(n-1)   Sn = a0(ratio^n -1)/(ratio -1)
    // -> Sn=an/ratio^(n-1)*(ratio^n -1)/(ratio -1), assume y=Sn/an
    // -> y*ratio^(n-1)=(ratio^n -1)/(ratio -1)
    // -> y*ratio^n - y*ratio^(n-1) -ratio^n +1 = 0

    // Solve using Newton-Raphson iteration to determine the expansion
    // ratio giving the desired overall timestep change.
    // Note max iteration count to avoid hanging; function generally
    // converges in 5 iterations or so.
    for (label iter = 0; iter < 100; iter++)
    {
        // Dimensionless equation
        scalar f = (y-1)*pow(ratio, n)+1-y*pow(ratio, n-1);
        scalar dfdratio = (y-1)*n*pow(ratio, n-1)-y*(n-1)*pow(ratio, n-2);
        scalar dratio = f/(dfdratio + SMALL);
        ratio -= dratio;
        // Exit if function is satisfied
        if (mag(f) < 1e-6)
        {
            return ratio;
        }
    }

    if (debug)
    {
        WarningInFunction
            << "Did not converge to find new timestep growth factor given "
            << "overall factor " << y << " and " << n << " steps to do it in."
            << endl << "    Returning current best guess " << ratio << endl;
    }
    return ratio;
}


void Foam::functionObjects::timeControl::calcDeltaTCoeff
(
    scalar& requiredDeltaTCoeff,
    const scalar wantedDT,
    const label nSteps,
    const scalar presentTime,
    const scalar timeToNextWrite,
    const bool rampDirectionUp
)
{
    const scalar writeInterval = writeControl_.interval();


    // Set new deltaT approximated to last step value
    // adjust to include geometric series sum of nSteps terms
    scalar newDeltaT = deltaT0_;

    if (seriesDTCoeff_ != GREAT)
    {
        newDeltaT *= seriesDTCoeff_;
    }

    // Recalculate new deltaTCoeff_ based on rounded steps
    requiredDeltaTCoeff = Foam::exp(Foam::log(wantedDT/newDeltaT)/nSteps);

    // Calculate total time required with given dT increment
    // to reach specified dT value
    scalar requiredTimeInterval = newDeltaT;

    // For nSteps = 1 during we get coeff = 1.0. Adding SMALL in denominator
    // might lead to SMALL requiredTimeInterval, and further unnecessary
    // calculations
    if (requiredDeltaTCoeff != 1.0)
    {
        requiredTimeInterval *=
            (pow(requiredDeltaTCoeff, nSteps) - 1)
           /(requiredDeltaTCoeff - 1);
    }

    // Nearest time which is multiple of wantedDT, after
    // requiredTimeInterval
    scalar timeToNextMultiple = -presentTime;

    if (rampDirectionUp)
    {
        timeToNextMultiple +=
            ceil
            (
                (presentTime + requiredTimeInterval)
               /wantedDT
            )
           *wantedDT;
    }
    else
    {
        timeToNextMultiple +=
            floor
            (
                (presentTime - requiredTimeInterval)
               /wantedDT
            )
           *wantedDT;
    }

    // Check nearest time and decide parameters
    if (timeToNextWrite <= timeToNextMultiple)
    {
        // As the ramping can't be achieved, use current, unadapted deltaT
        // (from e.g. maxCo calculation in solver) and adjust when
        // writeInterval is closer than this deltaT
        newDeltaT = deltaT0_;
        if (timeToNextWrite < newDeltaT)
        {
            newDeltaT = timeToNextWrite;
        }

        // Corner case when required time span is less than writeInterval
        if (requiredTimeInterval > writeInterval)
        {
            WarningInFunction
                << "With given ratio needed time span "
                << requiredTimeInterval
                << " exceeds available writeInterval "
                << writeInterval << nl
                << "Disabling all future time step ramping"
                << endl;
            deltaTCoeff_ = GREAT;
            newDeltaT = wantedDT;
        }
        else
        {
            // Reset the series ratio, as this could change later
            seriesDTCoeff_ = GREAT;
            if (debug)
            {
                Info<< "Disabling ramping until time "
                    << presentTime+timeToNextWrite << endl;
            }
        }
        requiredDeltaTCoeff = newDeltaT/deltaT0_;
    }
    else
    {
        // Find the minimum number of time steps (with constant deltaT
        // variation) to get us to the writeInterval and multiples of wantedDT.
        // Increases the number of nSteps to do the adjustment in until it
        // has reached one that is possible. Returns the new expansionratio.

        scalar y = timeToNextMultiple/wantedDT;
        label requiredSteps = nSteps;

        scalar ratioEstimate = deltaTCoeff_;
        scalar ratioMax = deltaTCoeff_;

        if (seriesDTCoeff_ != GREAT)
        {
            ratioEstimate = seriesDTCoeff_;
        }

        if (!rampDirectionUp)
        {
            ratioEstimate = 1/ratioEstimate;
            ratioMax = 1/ratioMax;
            requiredSteps *= -1;
        }
        // Provision for fall-back value if we can't satisfy requirements
        bool searchConverged = false;
        for (label iter = 0; iter < 100; iter++)
        {
            // Calculate the new expansion and the overall time step changes
            const scalar newRatio = calcExpansion
            (
                ratioEstimate,
                y,
                requiredSteps
            );
            const scalar deltaTLoop =
                wantedDT
              / pow(newRatio, mag(requiredSteps)-1);
            scalar firstDeltaRatio = deltaTLoop/deltaT0_;
            // Avoid division by zero for ratio = 1.0
            scalar Sn =
                deltaTLoop
               *(pow(newRatio, mag(requiredSteps))-1)
               /(newRatio-1+SMALL);

            if (debug)
            {
                Info<< " nSteps " << requiredSteps
                    << " ratioEstimate " << ratioEstimate
                    << " a0 " << deltaTLoop
                    << " firstDeltaRatio " << firstDeltaRatio
                    << " Sn " << Sn << " Sn+Time " << Sn+presentTime
                    << " seriesRatio "<< newRatio << endl;
            }

            // If series ratio becomes unity check next deltaT multiple
            // Do the same if first ratio is decreasing when ramping up!
            if
            (
                (firstDeltaRatio < 1.0 && rampDirectionUp)
             || (firstDeltaRatio > 1.0 && !rampDirectionUp)
             || newRatio == 1.0
            )
            {
                y += 1;
                requiredSteps = mag(nSteps);
                if (debug)
                {
                    Info<< "firstDeltaRatio " << firstDeltaRatio << " rampDir"
                        << rampDirectionUp << " newRatio " << newRatio
                        << " y " << y << " steps " << requiredSteps <<endl;
                }
                continue;
            }

            if (firstDeltaRatio > ratioMax && rampDirectionUp)
            {
                requiredSteps++;
                if (debug)
                {
                    Info<< "First ratio " << firstDeltaRatio
                        << " exceeds threshold " << ratioMax << nl
                        << "Increasing required steps " << requiredSteps
                        << endl;
                }
            }
            else if (firstDeltaRatio < ratioMax && !rampDirectionUp)
            {
                requiredSteps--;
                if (debug)
                {
                    Info<< "First ratio " << firstDeltaRatio
                        << " exceeds threshold " << ratioMax << nl
                        << "Decreasing required steps " << requiredSteps
                        << endl;
                }
            }
            else if
            (
                 (newRatio > ratioMax && rampDirectionUp)
             ||  (newRatio < ratioMax && !rampDirectionUp)
            )
            {
                y += 1;
                requiredSteps = nSteps;
                if (debug)
                {
                    Info<< "Series ratio " << newRatio << "exceeds threshold "
                        << ratioMax << nl << "Consider next deltaT multiple "
                        << y*wantedDT + presentTime << endl;
                }
            }
            else
            {
                seriesDTCoeff_ = newRatio;
                // Reset future series expansion at last step
                if (requiredSteps == 1)
                {
                    Sn = deltaT0_*firstDeltaRatio;
                    seriesDTCoeff_ = GREAT;
                }
                // Situations where we achieve wantedDT but fail to achieve
                // multiple of writeInterval
                scalar jumpError =
                   remainder(Sn + presentTime, wantedDT) / wantedDT;

                if (mag(jumpError) > ROOTSMALL)
                {
                    requiredSteps = label(timeToNextWrite/wantedDT);
                    firstDeltaRatio = timeToNextWrite/requiredSteps/deltaT0_;
                }
                if (debug)
                {
                    Info<< "All conditions satisfied deltaT0_:" << deltaT0_
                        << " calculated deltaTCoeff:" << firstDeltaRatio
                        << " series ratio set to:" << seriesDTCoeff_ << endl;
                }
                searchConverged = true;
                requiredDeltaTCoeff = firstDeltaRatio;
                break;
            }
        }

        if (!searchConverged)
        {
            if (debug)
            {
                // Not converged. Return fall-back value
                WarningInFunction
                    << "Did not converge to find new timestep growth factor"
                    << " given overall factor " << y
                    << " and " << requiredSteps << " steps to do it in."
                    << nl << "    Falling back to non-adjusted deltatT "
                    << deltaT0_ << endl;
            }
            requiredDeltaTCoeff = 1.0;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::timeControl::timeControl
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    dict_(dict),
    controlMode_(controlMode::TIME),
    timeStart_(-VGREAT),
    timeEnd_(VGREAT),
    triggerStart_(labelMin),
    triggerEnd_(labelMax),
    nStepsToStartTimeChange_(labelMax),
    executeControl_(runTime, dict, "execute"),
    writeControl_(runTime, dict, "write"),
    foPtr_(functionObject::New(name, runTime, dict_)),
    executeTimeIndex_(-1),
    deltaT0_(0),
    seriesDTCoeff_(GREAT)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeControl::entriesPresent(const dictionary& dict)
{
    if
    (
        Foam::timeControl::entriesPresent(dict, "write")
     || Foam::timeControl::entriesPresent(dict, "output") // backwards compat
     || Foam::timeControl::entriesPresent(dict, "execute")
     || dict.found("timeStart")
     || dict.found("timeEnd")
     || dict.found("triggerStart")
     || dict.found("triggerEnd")
    )
    {
        return true;
    }

    return false;
}


bool Foam::functionObjects::timeControl::execute()
{
    // Store old deltaT for reference
    deltaT0_ = time_.deltaTValue();

    if (active() && (postProcess || executeControl_.execute()))
    {
        executeTimeIndex_ = time_.timeIndex();
        foPtr_->execute();
    }

    return true;
}


bool Foam::functionObjects::timeControl::execute(const label subIndex)
{
    if (active())
    {
        // Call underlying function object directly
        foPtr_->execute(subIndex);
    }

    return true;
}


bool Foam::functionObjects::timeControl::write()
{
    if (active() && (postProcess || writeControl_.execute()))
    {
        // Ensure written results reflect the current state
        if (executeTimeIndex_ != time_.timeIndex())
        {
            executeTimeIndex_ = time_.timeIndex();
            foPtr_->execute();
        }

        foPtr_->write();
    }

    return true;
}


bool Foam::functionObjects::timeControl::end()
{
    if (active() && (executeControl_.execute() || writeControl_.execute()))
    {
        foPtr_->end();
    }

    return true;
}


bool Foam::functionObjects::timeControl::adjustTimeStep()
{
    if (active())
    {
        if
        (
            writeControl_.control()
         == Foam::timeControl::ocAdjustableRunTime
        )
        {
            const scalar presentTime = time_.value();

            //const scalar oldDeltaT = time_.deltaTValue();

            // Call underlying fo to do time step adjusting
            foPtr_->adjustTimeStep();
            const scalar wantedDT = time_.deltaTValue();

            //Pout<< "adjustTimeStep by " << foPtr_->name()
            //    << " from " << oldDeltaT << " to " << wantedDT << endl;

            const label writeTimeIndex = writeControl_.executionIndex();
            const scalar writeInterval = writeControl_.interval();

            // Get time to next writeInterval
            scalar timeToNextWrite =
                (writeTimeIndex+1)*writeInterval
              - (presentTime - time_.startTime().value());
            if (timeToNextWrite <= 0.0)
            {
                timeToNextWrite = writeInterval;
            }

            // The target DT (coming from fixed dT or maxCo etc) may not be
            // an integral multiple of writeInterval. If so ignore any step
            // calculation and give priority to writeInterval.
            // Note looser tolerance on the relative intervalError - SMALL is
            // too strict.
            scalar intervalError =
                remainder(writeInterval, wantedDT)/writeInterval;
            if
            (
                 mag(intervalError) > ROOTSMALL
              || deltaTCoeff_ == GREAT
            )
            {
                scalar deltaT = time_.deltaTValue();
                scalar nSteps = timeToNextWrite/deltaT;

                // For tiny deltaT the label can overflow!
                if
                (
                    nSteps < scalar(labelMax)
                 && (
                        deltaTCoeff_ != GREAT
                     || nSteps < nStepsToStartTimeChange_
                    )
                )
                {
                    // nSteps can be < 1 so make sure at least 1

                    label nStepsToNextWrite = max(1, round(nSteps));
                    scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

                    // Backwards compatibility: clip deltaT to 0.2 .. 2
                    scalar clipThreshold = 2;
                    if (deltaTCoeff_ != GREAT)
                    {
                        clipThreshold = deltaTCoeff_;
                    }
                    // Adjust time step
                    if (newDeltaT >= deltaT)
                    {
                        deltaT = min(newDeltaT, clipThreshold*deltaT);
                    }
                    else
                    {
                        clipThreshold = 1/clipThreshold;
                        deltaT = max(newDeltaT, clipThreshold*deltaT);
                    }

                    const_cast<Time&>(time_).setDeltaT(deltaT, false);
                }
            }
            else
            {
                // Set initial approximation to user defined ratio.
                scalar requiredDeltaTCoeff = deltaTCoeff_;

                // Use if we have series expansion ratio from last attempt.
                // Starting from user defined coefficient, might get different
                // steps and convergence time (could be recursive - worst case)
                // This is more likely to happen when coefficient limit is high
                if (seriesDTCoeff_ != GREAT)
                {
                    requiredDeltaTCoeff = seriesDTCoeff_;
                }
                // Avoid divide by zero if we need ratio = 1
                if (1/Foam::log(requiredDeltaTCoeff) > scalar(labelMax))
                {
                    requiredDeltaTCoeff = deltaTCoeff_;
                }

                // Try and observe the deltaTCoeff and the writeInterval
                // and the wanted deltaT. Below code calculates the
                // minimum number of time steps in which to end up on
                // writeInterval+n*deltaT with the minimum possible deltaTCoeff.
                // This is non-trivial!

                // Approx steps required from present dT to required dT
                label nSteps = 0;
                if (wantedDT < deltaT0_)
                {
                    nSteps = label
                    (
                        floor
                        (
                            Foam::log(wantedDT/deltaT0_)
                           /Foam::log(requiredDeltaTCoeff + SMALL)
                        )
                    );
                }
                else
                {
                    nSteps = label
                    (
                        ceil
                        (
                            Foam::log(wantedDT/deltaT0_)
                           /Foam::log(requiredDeltaTCoeff + SMALL)
                        )
                    );
                }
                // Negative steps -> ramp down needed,
                // zero steps -> we are at writeInterval-multiple time!
                if (nSteps < 1)
                {
                    // Not enough steps to do gradual time step adjustment
                    // if earlier deltaT larger than wantedDT
                    if (wantedDT < deltaT0_)
                    {
                        requiredDeltaTCoeff = 1/requiredDeltaTCoeff;
                        calcDeltaTCoeff
                        (
                            requiredDeltaTCoeff,
                            wantedDT,
                            nSteps,
                            presentTime,
                            timeToNextWrite,
                            false
                        );
                    }
                    else
                    {
                        if (timeToNextWrite > wantedDT)
                        {
                            requiredDeltaTCoeff = wantedDT/deltaT0_;
                        }
                        else
                        {
                            requiredDeltaTCoeff = timeToNextWrite/deltaT0_;
                        }
                    }
                    // When allowed deltaTCoeff can't give required ramp
                    if (requiredDeltaTCoeff > deltaTCoeff_ && debug)
                    {
                        WarningInFunction
                            << "Required deltaTCoeff "<< requiredDeltaTCoeff
                            << " is larger than allowed value "<< deltaTCoeff_
                            << endl;
                    }
                }
                else
                {
                    calcDeltaTCoeff
                    (
                        requiredDeltaTCoeff,
                        wantedDT,
                        nSteps,
                        presentTime,
                        timeToNextWrite,
                        true
                    );
                }

                // Calculate deltaT to be used based on ramp
                scalar newDeltaT = deltaT0_*requiredDeltaTCoeff;

                // Set time, deltaT adjustments for writeInterval purposes
                // are already taken care. Hence disable auto-update
                const_cast<Time&>(time_).setDeltaT(newDeltaT, false);

                if (seriesDTCoeff_ < 1.0)
                {
                    requiredDeltaTCoeff = 1/requiredDeltaTCoeff;
                    seriesDTCoeff_ = 1/seriesDTCoeff_;
                }
            }
        }
        else
        {
            foPtr_->adjustTimeStep();
            const scalar wantedDT = time_.deltaTValue();

            if (deltaTCoeff_ != GREAT)
            {
                // Clip time step change to deltaTCoeff

                scalar requiredDeltaTCoeff =
                (
                    max
                    (
                        1.0/deltaTCoeff_,
                        min(deltaTCoeff_, wantedDT/deltaT0_)
                    )
                );

                scalar newDeltaT = deltaT0_*requiredDeltaTCoeff;

                // Set time, deltaT adjustments for writeInterval purposes
                // are already taken care. Hence disable auto-update
                const_cast<Time&>( time_).setDeltaT(newDeltaT, false);
            }
            else
            {
                // Set the DT without any ramping
                const_cast<Time&>( time_).setDeltaT(wantedDT, false);
            }
        }
    }

    return true;
}


bool Foam::functionObjects::timeControl::read(const dictionary& dict)
{
    if (dict != dict_)
    {
        dict_ = dict;

        writeControl_.read(dict);
        executeControl_.read(dict);
        readControls();

        // Forward to underlying function object
        return foPtr_->read(dict);
    }

    return false;
}


bool Foam::functionObjects::timeControl::filesModified() const
{
    return active() ? foPtr_->filesModified() : false;
}


void Foam::functionObjects::timeControl::updateMesh(const mapPolyMesh& mpm)
{
    if (active())
    {
        foPtr_->updateMesh(mpm);
    }
}


void Foam::functionObjects::timeControl::movePoints(const polyMesh& mesh)
{
    if (active())
    {
        foPtr_->movePoints(mesh);
    }
}


// ************************************************************************* //
