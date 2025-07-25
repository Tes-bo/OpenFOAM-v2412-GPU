/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::functionObjects::volRegion::useAllCells() const noexcept
{
    return (regionType_ == vrtAll);
}


inline Foam::functionObjects::volRegion::regionTypes
Foam::functionObjects::volRegion::regionType() const noexcept
{
    return regionType_;
}


inline Foam::label Foam::functionObjects::volRegion::nCells() const
{
    #ifdef FULLDEBUG
    if (requireUpdate_)
    {
        FatalErrorInFunction
            << "Retrieving cached values that are not up-to-date" << nl
            << exit(FatalError);
    }
    #endif

    return nCells_;
}


inline Foam::scalar Foam::functionObjects::volRegion::V() const
{
    #ifdef FULLDEBUG
    if (requireUpdate_)
    {
        FatalErrorInFunction
            << "Retrieving cached values that are not up-to-date" << nl
            << exit(FatalError);
    }
    #endif

    return V_;
}


// ************************************************************************* //
