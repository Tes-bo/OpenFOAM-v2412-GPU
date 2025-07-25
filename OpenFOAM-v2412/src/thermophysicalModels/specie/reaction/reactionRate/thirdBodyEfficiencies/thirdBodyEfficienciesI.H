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

#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const scalarList& efficiencies
)
:
    scalarList(efficiencies),
    species_(species)
{
    if (size() != species_.size())
    {
        FatalErrorInFunction
            << "Number of efficiencies = " << size()
            << " is not equal to the number of species " << species_.size()
            << exit(FatalError);
    }
}


inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const dictionary& dict
)
:
    scalarList(species.size()),
    species_(species)
{
    List<Tuple2<word, scalar>> coeffs;
    if (dict.readIfPresent("coeffs", coeffs))
    {
        if (coeffs.size() != species_.size())
        {
            FatalIOErrorInFunction(dict)
                << "Number of efficiencies = " << coeffs.size()
                << " is not equal to the number of species " << species_.size()
                << exit(FatalIOError);
        }

        for (const auto& coeff : coeffs)
        {
            const label idx = species.find(coeff.first());

            if (idx < 0)
            {
                FatalIOErrorInFunction(dict)
                    << "No such species: " << coeffs.first()
                    << ", available: " << flatOutput(species) << endl
                    << exit(FatalIOError);
            }

            operator[](idx) = coeff.second();
        }
    }
    else
    {
        const scalar defaultEff = dict.get<scalar>("defaultEfficiency");
        scalarList::operator=(defaultEff);
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyEfficiencies::M(const scalarList& c) const
{
    scalar M = 0.0;
    forAll(*this, i)
    {
        M += operator[](i)*c[i];
    }

    return M;
}


inline void Foam::thirdBodyEfficiencies::write(Ostream& os) const
{
    List<Tuple2<word, scalar>> coeffs(species_.size());
    forAll(coeffs, i)
    {
        coeffs[i].first() = species_[i];
        coeffs[i].second() = operator[](i);
    }

    os.writeEntry("coeffs", coeffs);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thirdBodyEfficiencies& tbes
)
{
    tbes.write(os);
    return os;
}


// ************************************************************************* //
