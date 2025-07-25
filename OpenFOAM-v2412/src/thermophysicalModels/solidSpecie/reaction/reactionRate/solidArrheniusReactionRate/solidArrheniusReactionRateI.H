/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::solidArrheniusReactionRate::solidArrheniusReactionRate
(
    const scalar A,
    const scalar Ta,
    const scalar Tcrit
)
:
    A_(A),
    Ta_(Ta),
    Tcrit_(Tcrit)
{}


inline Foam::solidArrheniusReactionRate::solidArrheniusReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.get<scalar>("A")),
    Ta_(dict.get<scalar>("Ta")),
    Tcrit_(dict.get<scalar>("Tcrit"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::solidArrheniusReactionRate::operator()
(
    const scalar,
    const scalar T,
    const scalarField&
) const
{
    if (T < Tcrit_)
    {
        return 0;
    }

    return A_*exp(-Ta_/T);
}


inline void Foam::solidArrheniusReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("Tcrit", Tcrit_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const solidArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
