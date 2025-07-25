/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "liquid.H"
#include "NoneFunction1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquid, 0);
    addToRunTimeSelectionTable(liquidProperties, liquid, dictionary);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static autoPtr<Function1<Type>> NewOrNone
(
    const word& entryName,
    const dictionary& dict
)
{
    autoPtr<Function1<Type>> ptr
    (
        Function1<Type>::NewIfPresent(entryName, dict)
    );

    if (!ptr)
    {
        ptr.reset
        (
            new Function1Types::None<Type>(entryName, dict)
        );
    }

    return ptr;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquid::liquid(const dictionary& dict)
:
    liquidProperties(dict),
    rho_(NewOrNone<scalar>("rho", dict)),
    pv_(NewOrNone<scalar>("pv", dict)),
    hl_(NewOrNone<scalar>("hl", dict)),
    Cp_(NewOrNone<scalar>("Cp", dict)),
    h_(NewOrNone<scalar>("h", dict)),
    Cpg_(NewOrNone<scalar>("Cpg", dict)),
    B_(NewOrNone<scalar>("B", dict)),
    mu_(NewOrNone<scalar>("mu", dict)),
    mug_(NewOrNone<scalar>("mug", dict)),
    kappa_(NewOrNone<scalar>("kappa", dict)),
    kappag_(NewOrNone<scalar>("kappag", dict)),
    sigma_(NewOrNone<scalar>("sigma", dict)),
    D_(NewOrNone<scalar>("D", dict))
{}



Foam::liquid::liquid(const liquid& rhs)
:
    liquidProperties(rhs),
    rho_(rhs.rho_.clone()),
    pv_(rhs.pv_.clone()),
    hl_(rhs.hl_.clone()),
    Cp_(rhs.Cp_.clone()),
    h_(rhs.h_.clone()),
    Cpg_(rhs.Cpg_.clone()),
    B_(rhs.B_.clone()),
    mu_(rhs.mu_.clone()),
    mug_(rhs.mug_.clone()),
    kappa_(rhs.kappa_.clone()),
    kappag_(rhs.kappag_.clone()),
    sigma_(rhs.sigma_.clone()),
    D_(rhs.D_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquid::writeData(Ostream& os) const
{
    liquidProperties::writeData(os); os << nl;
    rho_->writeData(os); os << nl;
    pv_->writeData(os); os << nl;
    hl_->writeData(os); os << nl;
    Cp_->writeData(os); os << nl;
    h_->writeData(os); os << nl;
    Cpg_->writeData(os); os << nl;
    B_->writeData(os); os << nl;
    mu_->writeData(os); os << nl;
    mug_->writeData(os); os << nl;
    kappa_->writeData(os); os << nl;
    kappag_->writeData(os); os << nl;
    sigma_->writeData(os); os << nl;
    D_->writeData(os); os << endl;
}


// ************************************************************************* //
