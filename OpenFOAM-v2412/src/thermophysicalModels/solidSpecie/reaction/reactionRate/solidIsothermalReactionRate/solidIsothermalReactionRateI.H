/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::solidIsothermalReactionRate::solidIsothermalReactionRate
(
    const scalar C,
    const scalar Cp,
    const scalar Tpc,
    const scalar Elat
)
:
    C_(C),
    Cp_(Cp),
    Tpc_(Tpc),
    Elat_(Elat)
{}


inline Foam::solidIsothermalReactionRate::solidIsothermalReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    C_(dict.get<scalar>("C")),
    Cp_(dict.get<scalar>("Cp")),
    Tpc_(dict.get<scalar>("Tpc")),
    Elat_(dict.get<scalar>("Elat"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidIsothermalReactionRate::operator()
(
    const scalar,
    const scalar T,
    const scalarField&
) const
{
    if (T < Tpc_)
    {
        return 0;
    }

    return -C_*Cp_*(T - Tpc_)/Elat_;
}


inline void Foam::solidIsothermalReactionRate::write(Ostream& os) const
{
    os.writeEntry("C", C_);
    os.writeEntry("Cp", Cp_);
    os.writeEntry("Tpc", Tpc_);
    os.writeEntry("Elat", Elat_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solidIsothermalReactionRate& iso
)
{
    iso.write(os);
    return os;
}


// ************************************************************************* //
