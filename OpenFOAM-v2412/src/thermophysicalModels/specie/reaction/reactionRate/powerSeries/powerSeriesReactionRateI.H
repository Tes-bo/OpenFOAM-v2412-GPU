/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

inline Foam::powerSeriesReactionRate::powerSeriesReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const FixedList<scalar, nCoeff_> coeffs
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta),
    coeffs_(coeffs)
{}


inline Foam::powerSeriesReactionRate::powerSeriesReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.get<scalar>("A")),
    beta_(dict.get<scalar>("beta")),
    Ta_(dict.get<scalar>("Ta")),
    coeffs_(dict.lookup("coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::powerSeriesReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&
) const
{
    scalar lta = A_;

    if (mag(beta_) > VSMALL)
    {
        lta *= pow(T, beta_);
    }

    scalar expArg = 0.0;

    forAll(coeffs_, n)
    {
        expArg += coeffs_[n]/pow(T, n + 1);
    }

    lta *= exp(expArg);

    return lta;
}


inline void Foam::powerSeriesReactionRate::write(Ostream& os) const
{
    os.writeEntry("A", A_);
    os.writeEntry("beta", beta_);
    os.writeEntry("Ta", Ta_);
    os.writeEntry("coeffs", coeffs_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const powerSeriesReactionRate& psrr
)
{
    psrr.write(os);
    return os;
}


// ************************************************************************* //
