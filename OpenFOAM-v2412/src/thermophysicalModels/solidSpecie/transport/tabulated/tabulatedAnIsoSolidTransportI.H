/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulatedAnIsoSolidTransport<Thermo>::tabulatedAnIsoSolidTransport
(
    const Thermo& t,
    const autoPtr<Function1<vector>>& kappa

)
:
    Thermo(t),
    kappa_(kappa.clone())
{}


template<class Thermo>
inline Foam::tabulatedAnIsoSolidTransport<Thermo>::tabulatedAnIsoSolidTransport
(
    const word& name,
    const tabulatedAnIsoSolidTransport& pt
)
:
    Thermo(name, pt),
    kappa_(pt.kappa_.clone())
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedAnIsoSolidTransport<Thermo>>
Foam::tabulatedAnIsoSolidTransport<Thermo>::clone() const
{
    return autoPtr<tabulatedAnIsoSolidTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulatedAnIsoSolidTransport<Thermo>>
Foam::tabulatedAnIsoSolidTransport<Thermo>::New(const dictionary& dict)
{
    return autoPtr<tabulatedAnIsoSolidTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulatedAnIsoSolidTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return 0;
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedAnIsoSolidTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return mag(kappa_->value(T));
}


template<class Thermo>
inline Foam::vector Foam::tabulatedAnIsoSolidTransport<Thermo>::Kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_->value(T);
}


template<class Thermo>
inline Foam::scalar Foam::tabulatedAnIsoSolidTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return kappa(p, T)/this->Cp(p, T);
}


// ************************************************************************* //
