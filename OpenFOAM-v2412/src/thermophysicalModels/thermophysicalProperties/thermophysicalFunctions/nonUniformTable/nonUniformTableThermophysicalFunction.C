/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "nonUniformTableThermophysicalFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonUniformTable, 0);

    addToRunTimeSelectionTable
    (
        thermophysicalFunction,
        nonUniformTable,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonUniformTable::nonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    values_(),
    Trange_(),
    deltaT_(GREAT),
    jumpTable_()
{
    dict.readEntry(name_, values_);

    if (values_.size() < 2)
    {
        FatalIOErrorInFunction(dict)
            << "Table" << nl
            << "    " << name_ << nl
            << "    has fewer than 2 entries." << nl
            << exit(FatalIOError);
    }

    Trange_.min() = values_.first().first();
    Trange_.max() = values_.last().first();

    for (label i = 1; i < values_.size(); ++i)
    {
        #ifdef FULLDEBUG
        // Check list is monotonically increasing...
        if (values_[i].first() <= values_[i-1].first())
        {
            FatalErrorInFunction
                << "Table" << nl
                << "    " << name_ << nl
                << "    out-of-order value: " << values_[i].first()
                << " at index " << i << nl
                << exit(FatalError);
        }
        #endif

        deltaT_ = min(deltaT_, values_[i].first() - values_[i-1].first());
    }

    deltaT_ *= 0.9;

    jumpTable_.resize(Trange_.mag()/deltaT_ + 1);

    label i = 0;
    forAll(jumpTable_, j)
    {
        const scalar T = Trange_.min() + j*deltaT_;

        if (T > values_[i+1].first())
        {
            ++i;
        }

        jumpTable_[j] = i;
    }
}


Foam::nonUniformTable::nonUniformTable
(
    const dictionary& dict
)
:
    nonUniformTable("values", dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::nonUniformTable::f
(
    scalar p,
    scalar T
) const
{
    const label i = index(p, T);
    const scalar Ti = values_[i].first();
    const scalar lambda = (T - Ti)/(values_[i + 1].first() - Ti);

    return
        values_[i].second()
      + lambda*(values_[i + 1].second() - values_[i].second());
}


Foam::scalar Foam::nonUniformTable::dfdT
(
    scalar p,
    scalar T
) const
{
    const label i = index(p, T);

    return
        (values_[i + 1].second() - values_[i].second())
       /(values_[i + 1].first() - values_[i].first());
}


void Foam::nonUniformTable::writeData(Ostream& os) const
{
    os.writeEntry("values", values_);
}


// ************************************************************************* //
