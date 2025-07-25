/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "WallSiteData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::WallSiteData<Type>::WallSiteData()
:
    patchi_(),
    wallData_()
{}


template<class Type>
Foam::WallSiteData<Type>::WallSiteData
(
    label patchi,
    const Type& wallData
)
:
    patchi_(patchi),
    wallData_(wallData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::WallSiteData<Type>::~WallSiteData()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::WallSiteData<Type>::operator==
(
    const WallSiteData<Type>& rhs
) const
{
    return patchi_ == rhs.patchi_ && wallData_ == rhs.wallData_;
}


template<class Type>
bool Foam::WallSiteData<Type>::operator!=
(
    const WallSiteData<Type>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    WallSiteData<Type>& wIS
)
{
    is  >> wIS.patchi_ >> wIS.wallData_;

    is.check(FUNCTION_NAME);
    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WallSiteData<Type>& wIS
)
{
    os  << wIS.patchi_ << token::SPACE << wIS.wallData_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
