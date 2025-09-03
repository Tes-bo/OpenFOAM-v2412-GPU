/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionObjects::fieldAverageItem::EXT_MEAN
(
    "Mean"
);


const Foam::word Foam::functionObjects::fieldAverageItem::EXT_PRIME2MEAN
(
    "Prime2Mean"
);


const Foam::Enum
<
    Foam::functionObjects::fieldAverageItem::baseType
>
Foam::functionObjects::fieldAverageItem::baseTypeNames_
({
    { baseType::ITER, "iteration" },
    { baseType::TIME, "time" },
});


const Foam::Enum
<
    Foam::functionObjects::fieldAverageItem::windowType
>
Foam::functionObjects::fieldAverageItem::windowTypeNames_
({
    { windowType::NONE, "none" },
    { windowType::APPROXIMATE, "approximate" },
    { windowType::EXACT, "exact" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverageItem::fieldAverageItem()
:
    active_(false),
    fieldName_("unknown"),
    mean_(false),
    meanFieldName_("unknown"),
    prime2Mean_(false),
    prime2MeanFieldName_("unknown"),
    base_(baseType::ITER),
    totalIter_(0),
    totalTime_(-1),
    window_(-1),
    windowName_(""),
    windowType_(windowType::NONE),

    windowTimes_(),
    windowFieldNames_(),
    allowRestart_(true)
{}


Foam::functionObjects::fieldAverageItem::fieldAverageItem
(
    const fieldAverageItem& faItem
)
:
    active_(faItem.active_),
    fieldName_(faItem.fieldName_),
    mean_(faItem.mean_),
    meanFieldName_(faItem.meanFieldName_),
    prime2Mean_(faItem.prime2Mean_),
    prime2MeanFieldName_(faItem.prime2MeanFieldName_),
    base_(faItem.base_),
    totalIter_(faItem.totalIter_),
    totalTime_(faItem.totalTime_),
    window_(faItem.window_),
    windowName_(faItem.windowName_),
    windowType_(faItem.windowType_),

    windowTimes_(faItem.windowTimes_),
    windowFieldNames_(faItem.windowFieldNames_),
    allowRestart_(faItem.allowRestart_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverageItem::~fieldAverageItem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldAverageItem::addToWindow
(
    const word& fieldName,
    const scalar deltaT
)
{
    windowTimes_.push(deltaT);
    windowFieldNames_.push(fieldName);
}


void Foam::functionObjects::fieldAverageItem::evolve(const objectRegistry& obr)
{
    totalIter_++;
    totalTime_ += obr.time().deltaTValue();
    forAllIters(windowTimes_, timeIter)
    {
        timeIter() += obr.time().deltaTValue();
    }

    // Remove any fields that have passed out of the window
    bool removeItem = true;

    while (removeItem && windowTimes_.size())
    {
        removeItem = !(inWindow(windowTimes_.first()));

        if (removeItem)
        {
            windowTimes_.pop();
            const word fieldName = windowFieldNames_.pop();

            //Info<< "evolve: removing field: " << fieldName << endl;
            obr.checkOut(fieldName);
        }
    }
}


void Foam::functionObjects::fieldAverageItem::clear
(
    const objectRegistry& obr,
    bool fullClean
)
{
    if (mean_)
    {
        obr.checkOut(meanFieldName_);
    }

    if (prime2Mean_)
    {
        obr.checkOut(prime2MeanFieldName_);
    }

    for (const word& fieldName : windowFieldNames_)
    {
        obr.checkOut(fieldName);
    }

    if (totalTime_ < 0 || fullClean)
    {
        totalIter_ = 0;
        totalTime_ = 0;
        windowTimes_.clear();
        windowFieldNames_.clear();
    }
}


bool Foam::functionObjects::fieldAverageItem::readState(const dictionary& dict)
{
    dict.readEntry("totalIter", totalIter_);
    dict.readEntry("totalTime", totalTime_);

    if (window_ > 0)
    {
        dict.readEntry("windowTimes", windowTimes_);
        dict.readEntry("windowFieldNames", windowFieldNames_);
    }

    return true;
}


void Foam::functionObjects::fieldAverageItem::writeState
(
    dictionary& dict
) const
{
    dict.add("totalIter", totalIter_);
    dict.add("totalTime", totalTime_);

    if (window_ > 0)
    {
        dict.add("windowTimes", windowTimes_);
        dict.add("windowFieldNames", windowFieldNames_);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldAverageItem::operator=
(
    const fieldAverageItem& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    // Set updated values
    active_ = rhs.active_;
    fieldName_ = rhs.fieldName_;
    mean_ = rhs.mean_;
    meanFieldName_ = rhs.meanFieldName_;
    prime2Mean_ = rhs.prime2Mean_;
    prime2MeanFieldName_ = rhs.prime2MeanFieldName_;
    base_ = rhs.base_;
    totalIter_ = rhs.totalIter_;
    totalTime_ = rhs.totalTime_;
    window_ = rhs.window_;
    windowName_ = rhs.windowName_;
    windowType_ = rhs.windowType_;
    windowTimes_ = rhs.windowTimes_;
    windowFieldNames_ = rhs.windowFieldNames_;
    allowRestart_ = rhs.allowRestart_;
}


// ************************************************************************* //
