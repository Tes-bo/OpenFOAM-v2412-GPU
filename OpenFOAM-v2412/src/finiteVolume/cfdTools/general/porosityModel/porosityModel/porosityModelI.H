/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline const Foam::coordinateSystem& Foam::porosityModel::csys() const
{
    return *csysPtr_;
}


inline Foam::label Foam::porosityModel::fieldIndex(const label i) const
{
    return (csysPtr_->uniform() ? 0 : i);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::porosityModel::name() const
{
    return name_;
}


inline bool Foam::porosityModel::active() const
{
    return active_;
}


inline const Foam::dictionary& Foam::porosityModel::dict() const
{
    return dict_;
}


inline const Foam::labelList& Foam::porosityModel::cellZoneIDs() const
{
    return cellZoneIDs_;
}


// ************************************************************************* //
