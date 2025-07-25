/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "lineSearch.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lineSearch, 0);
    defineRunTimeSelectionTable(lineSearch, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::lineSearch::coeffsDict()
{
    return dict_.optionalSubDict(type() + "Coeffs");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineSearch::lineSearch
(
    const dictionary& dict,
    const Time& time,
    updateMethod& UpdateMethod
)
:
    dict_(dict),
    lineSearchDict_
    (
        IOobject
        (
            "lineSearch",
            time.timeName(),
            "uniform",
            time,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ),
    directionalDeriv_(Zero),
    direction_(0),
    oldMeritValue_(Zero),
    newMeritValue_(Zero),
    prevMeritDeriv_
    (
        lineSearchDict_.getOrDefault<scalar>("prevMeritDeriv", Zero)
    ),
    initialStep_(dict.getOrDefault<scalar>("initialStep", 1.)),
    minStep_(dict.getOrDefault<scalar>("minStep", 0.3)),
    step_(Zero),
    iter_(lineSearchDict_.getOrDefault<label>("iter", 0)),
    innerIter_(0),
    maxIters_(dict.getOrDefault<label>("maxIters", 4)),
    extrapolateInitialStep_
    (
        dict.getOrDefault<bool>
        (
            "extrapolateInitialStep",
            false
        )
    ),
    stepUpdate_(stepUpdate::New(dict)),
    updateMethod_(UpdateMethod)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lineSearch> Foam::lineSearch::New
(
    const dictionary& dict,
    const Time& time,
    updateMethod& UpdateMethod
)
{
    autoPtr<lineSearch> lineSrch(nullptr);

    const word modelType(dict.getOrDefault<word>("type", "none"));

    Info<< "lineSearch type : " << modelType << endl;

    if (modelType != "none")
    {
        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                "lineSearch",
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        lineSrch.reset((ctorPtr(dict, time, UpdateMethod)).ptr());
    }
    else
    {
        Info<< "No line search method specified. "
            << "Proceeding with constant step" << endl;
    }

    return lineSrch;
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::lineSearch::setDeriv(const scalar deriv)
{
    directionalDeriv_ = deriv;
    stepUpdate_->setDeriv(deriv);
}


void Foam::lineSearch::setNewDeriv(const scalar deriv)
{
    // Does nothing in base
}


void Foam::lineSearch::setNewMeritValue(const scalar value)
{
    newMeritValue_ = value;
    stepUpdate_->setNewMeritValue(value);
}


void Foam::lineSearch::setOldMeritValue(const scalar value)
{
    oldMeritValue_ = value;
    stepUpdate_->setOldMeritValue(value);
}


void Foam::lineSearch::reset()
{
    innerIter_ = 0;
    if (extrapolateInitialStep_ && iter_ != 0)
    {
        // step_ = 2*(oldMeritValue_-prevMeritValue_)/directionalDeriv_;
        // Interpolate in order to get same improvement with the previous
        // optimisation cycle
        step_ =
            clamp(step_*prevMeritDeriv_/directionalDeriv_, minStep_, scalar(1));
        Info<< "\n------- Computing initial step-------" << endl;
        Info<< "old dphi(0) "  << prevMeritDeriv_ << endl;
        Info<< "dphi(0) "      << directionalDeriv_ << endl;
        Info<< "Setting initial step value " << step_ << endl << endl;
    }
    else
    {
        step_ = initialStep_;
    }
}


void Foam::lineSearch::updateStep(const scalar newStep)
{
    step_ = newStep;
}


void Foam::lineSearch::updateCorrection(scalarField& correction)
{
    correction *= step_;
}


bool Foam::lineSearch::loop()
{
    const bool isRunning = innerIter_ < maxIters_;

    if (isRunning)
    {
        ++innerIter_;
    }

    return isRunning;
}


bool Foam::lineSearch::computeGradient() const
{
    return false;
}


void Foam::lineSearch::postUpdate()
{
    this->operator++();
}


Foam::lineSearch& Foam::lineSearch::operator++()
{
    ++iter_;
    prevMeritDeriv_ = directionalDeriv_;
    lineSearchDict_.add<scalar>("prevMeritDeriv", prevMeritDeriv_, true);
    lineSearchDict_.add<label>("iter", iter_, true);
    if (lineSearchDict_.time().writeTime())
    {
        lineSearchDict_.regIOobject::writeObject
        (
            IOstreamOption(IOstreamOption::ASCII),
            true
        );
    }

    return *this;
}


Foam::lineSearch& Foam::lineSearch::operator++(int)
{
    return operator++();
}


// ************************************************************************* //
