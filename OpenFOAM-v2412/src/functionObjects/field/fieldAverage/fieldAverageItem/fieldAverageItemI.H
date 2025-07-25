/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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


bool Foam::functionObjects::fieldAverageItem::active() const
{
    return active_;
}


bool& Foam::functionObjects::fieldAverageItem::active()
{
    return active_;
}


const Foam::word& Foam::functionObjects::fieldAverageItem::fieldName() const
{
    return fieldName_;
}


bool Foam::functionObjects::fieldAverageItem::mean() const
{
    return mean_;
}


bool& Foam::functionObjects::fieldAverageItem::mean()
{
    return mean_;
}


const Foam::word&
Foam::functionObjects::fieldAverageItem::meanFieldName() const
{
    return meanFieldName_;
}


void Foam::functionObjects::fieldAverageItem::setMeanFieldName
(
    const word& name
)
{
    meanFieldName_ = name;
}


bool Foam::functionObjects::fieldAverageItem::prime2Mean() const
{
    return prime2Mean_;
}


bool& Foam::functionObjects::fieldAverageItem::prime2Mean()
{
    return prime2Mean_;
}


const Foam::word&
Foam::functionObjects::fieldAverageItem::prime2MeanFieldName() const
{
    return prime2MeanFieldName_;
}


void Foam::functionObjects::fieldAverageItem::setPrime2MeanFieldName
(
    const word& name
)
{
    prime2MeanFieldName_ = name;
}


const Foam::word& Foam::functionObjects::fieldAverageItem::base() const
{
    return baseTypeNames_[base_];
}


Foam::label Foam::functionObjects::fieldAverageItem::totalIter() const
{
    return totalIter_;
}


Foam::scalar Foam::functionObjects::fieldAverageItem::totalTime() const
{
    return totalTime_;
}


Foam::scalar Foam::functionObjects::fieldAverageItem::window() const
{
    return window_;
}


const Foam::word& Foam::functionObjects::fieldAverageItem::windowName() const
{
    return windowName_;
}


const Foam::FIFOStack<Foam::scalar>&
Foam::functionObjects::fieldAverageItem::windowTimes() const
{
    return windowTimes_;
}


const Foam::FIFOStack<Foam::word>&
Foam::functionObjects::fieldAverageItem::windowFieldNames() const
{
    return windowFieldNames_;
}


bool Foam::functionObjects::fieldAverageItem::allowRestart() const
{
    return allowRestart_;
}


Foam::scalar Foam::functionObjects::fieldAverageItem::dt
(
    const scalar deltaT
) const
{
    switch (base_)
    {
        case baseType::TIME:
        {
            return deltaT;
        }
        case baseType::ITER:
        {
            return scalar(1);
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << baseTypeNames_[base_]
                << abort(FatalError);
        }
    }

    return 0.0;
}


Foam::scalar Foam::functionObjects::fieldAverageItem::Dt() const
{
    switch (base_)
    {
        case baseType::TIME:
        {
            return totalTime_;
        }
        case baseType::ITER:
        {
            return scalar(totalIter_);
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << baseTypeNames_[base_]
                << abort(FatalError);
        }
    }

    return 0.0;
}


Foam::word Foam::functionObjects::fieldAverageItem::windowFieldName
(
    const word& prefix
) const
{
    return prefix + ':' + fieldName_ + ':' + Foam::name(totalIter_);
}


bool Foam::functionObjects::fieldAverageItem::inWindow(const scalar t) const
{
    switch (base_)
    {
        case baseType::ITER:
        {
            return round(t) <= round(window_) + 1;
        }
        case baseType::TIME:
        {
            return t <= window_;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled baseType enumeration "
                << baseTypeNames_[base_]
                << abort(FatalError);
        }
    }

    return false;
}


bool Foam::functionObjects::fieldAverageItem::storeWindowFields() const
{
    return windowType_ == windowType::EXACT;
}


bool Foam::functionObjects::fieldAverageItem::writeWindowFields() const
{
    return (allowRestart_ && window_ > 0);
}


// ************************************************************************* //
