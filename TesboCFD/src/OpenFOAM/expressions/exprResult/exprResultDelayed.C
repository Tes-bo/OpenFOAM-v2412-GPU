/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
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

#include "exprResultDelayed.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "addToRunTimeSelectionTable.H"

// #include <cassert>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
    defineTypeName(exprResultDelayed);

    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultDelayed,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultDelayed,
        empty
    );

} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::exprResultDelayed::exprResultDelayed()
:
    exprResult(),
    name_("none"),
    startExpr_(),
    settingResult_(),
    storeInterval_(1),
    delay_(10)
{}


Foam::expressions::exprResultDelayed::exprResultDelayed
(
    const exprResultDelayed& rhs
)
:
    exprResult(rhs),
    name_(rhs.name_),
    startExpr_(rhs.startExpr_),
    settingResult_(rhs.settingResult_),
    storedValues_(rhs.storedValues_),
    storeInterval_(rhs.storeInterval_),
    delay_(rhs.delay_)
{}


Foam::expressions::exprResultDelayed::exprResultDelayed
(
    const dictionary& dict
)
:
    exprResult(dict.subOrEmptyDict("value")),
    name_(dict.get<word>("name")),
    startExpr_("startupValue", dict),
    storeInterval_(dict.get<scalar>("storeInterval")),
    delay_(dict.get<scalar>("delay"))
{
    const entry *eptr = dict.findEntry("storedValues", keyType::LITERAL);

    if (eptr && eptr->isStream())
    {
        storedValues_ = DLList<ValueAtTime>(eptr->stream());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::exprResultDelayed::updateReadValue
(
    const scalar& timeVal
)
{
    if (storedValues_.empty())
    {
        return false;
    }

    if (storedValues_.front().first() > (timeVal-delay_))
    {
        // No matching data yet
        return false;
    }

    if (storedValues_.size() <= 1)
    {
        FatalErrorInFunction
            << "Only one stored value at time " << timeVal
            << " for delayedVariable " << name() << nl
            << "Check the values for the interval " << storeInterval_
            << " and delay " << delay_ << nl
            << "Probably the interval is too large" << nl << endl
            << exit(FatalError);
    }

    auto current = storedValues_.cbegin();
    auto next = current;
    ++next;

    // The time without the delay offset
    const scalar newTime = (timeVal - delay_);

    while (next != storedValues_.end())
    {
        if (newTime >= current().first() && newTime <= next().first())
        {
            break;
        }

        current = next;
        ++next;
    }

    const scalar f =
    (
        (newTime - current().first())
      / (next().first() - current().first())
    );

    exprResult val((1-f)*current().second() + f*next().second());

    setReadValue(val);

    return true;
}


void Foam::expressions::exprResultDelayed::setReadValue
(
    const exprResult& val
)
{
    exprResult::operator=(val);
}


void Foam::expressions::exprResultDelayed::storeValue
(
    const scalar& currTime
)
{
    bool append = storedValues_.empty();

    if (!append)
    {
        const scalar lastTime = storedValues_.back().first();

        if (lastTime + SMALL >= currTime)
        {
            // Times are essentially identical - replace value
        }
        else if ((currTime - lastTime) >= 0.999*storeInterval_)
        {
            append = true;
        }
        else
        {
            // Cannot store in the middle - abandon the attempt
            return;
        }
    }

    if (append)
    {
        // Append value

        const scalar oldLastTime =
        (
            storedValues_.empty()
          ? 0
          : storedValues_.back().first()
        );

        storedValues_.push_back(ValueAtTime(currTime, settingResult_));

        while
        (
            storedValues_.size() > 1
         && (oldLastTime - storedValues_.front().first()) >= delay_
        )
        {
            // Remove values that are older than delay_
            storedValues_.pop_front();
        }
    }
    else
    {
        // Replace value

        storedValues_.back().second() = settingResult_;
    }
}


void Foam::expressions::exprResultDelayed::writeDict(Ostream& os) const
{
    os.beginBlock();

    os.writeEntry("name", name_);

    os.writeEntry("startupValue", startExpr_);

    if (!settingResult_.valueType().empty())
    {
        os.writeEntry("settingResult", settingResult_);
    }

    os.writeEntry("storedValues", storedValues_);
    os.writeEntry("storeInterval", storeInterval_);
    os.writeEntry("delay", delay_);

    os.writeKeyword("value");
    os << static_cast<const exprResult&>(*this);

    os.endBlock();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::expressions::exprResultDelayed::operator=
(
    const exprResultDelayed& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    exprResult::operator=(rhs);

    name_ = rhs.name_;
    startExpr_ = rhs.startExpr_;
    settingResult_ = rhs.settingResult_;
    storedValues_ = rhs.storedValues_;
    storeInterval_ = rhs.storeInterval_;
    delay_ = rhs.delay_;
}


void Foam::expressions::exprResultDelayed::operator=
(
    const exprResult& rhs
)
{
    settingResult_ = rhs;
}


void Foam::expressions::exprResultDelayed::operator=
(
    exprResult&& rhs
)
{
    settingResult_ = std::move(rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    expressions::exprResultDelayed& data
)
{
    dictionary dict(is);

    data = expressions::exprResultDelayed(dict);

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const expressions::exprResultDelayed& data
)
{
    data.writeDict(os);
    return os;
}


// ************************************************************************* //
