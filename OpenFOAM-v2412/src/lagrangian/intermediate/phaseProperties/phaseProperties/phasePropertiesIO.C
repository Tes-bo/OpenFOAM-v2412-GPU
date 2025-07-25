/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "phaseProperties.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseProperties::phaseProperties(Istream& is)
:
    phaseProperties()
{
    is >> *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, phaseProperties& pp)
{
    is.check(FUNCTION_NAME);

    const dictionaryEntry phaseInfo(dictionary::null, is);
    const dictionary& dict = phaseInfo.dict();

    pp.phase_ = pp.phaseTypeNames[phaseInfo.keyword()];
    pp.stateLabel_ = pp.phaseToStateLabel(pp.phase_);

    // The names in the dictionary order of appearance
    pp.names_ = dict.toc();

    const label nComponents = pp.names_.size();

    pp.Y_.resize(nComponents, Zero);
    pp.carrierIds_.resize(nComponents, -1);

    for (label cmpti = 0; cmpti < nComponents; ++cmpti)
    {
        pp.Y_[cmpti] = dict.get<scalar>(pp.names_[cmpti]);
    }
    pp.checkTotalMassFraction();

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const phaseProperties& pp)
{
    os.check(FUNCTION_NAME);

    os.beginBlock(pp.phaseTypeNames[pp.phase_]);

    forAll(pp.names_, cmpti)
    {
        os.writeEntry(pp.names_[cmpti], pp.Y_[cmpti]);
    }

    os.endBlock();

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
