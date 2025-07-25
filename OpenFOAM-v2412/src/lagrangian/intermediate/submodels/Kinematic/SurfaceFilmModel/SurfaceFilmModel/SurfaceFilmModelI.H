/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "SurfaceFilmModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const Foam::dimensionedVector& Foam::SurfaceFilmModel<CloudType>::g()
const noexcept
{
    return g_;
}


template<class CloudType>
Foam::label Foam::SurfaceFilmModel<CloudType>::nParcelsTransferred()
const noexcept
{
    return nParcelsTransferred_;
}

template<class CloudType>
Foam::label& Foam::SurfaceFilmModel<CloudType>::nParcelsTransferred()
noexcept
{
    return nParcelsTransferred_;
}


template<class CloudType>
Foam::label Foam::SurfaceFilmModel<CloudType>::nParcelsInjected()
const noexcept
{
    return nParcelsInjected_;
}


template<class CloudType>
Foam::label& Foam::SurfaceFilmModel<CloudType>::nParcelsInjected()
noexcept
{
    return nParcelsInjected_;
}


template<class CloudType>
Foam::scalar Foam::SurfaceFilmModel<CloudType>::totalMassTransferred()
const noexcept
{
    return totalMassTransferred_;
}


template<class CloudType>
Foam::scalar& Foam::SurfaceFilmModel<CloudType>::totalMassTransferred()
noexcept
{
    return totalMassTransferred_;
}


// ************************************************************************* //
