/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "RanzMarshall.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RanzMarshall<CloudType>::RanzMarshall
(
    const dictionary& dict,
    CloudType& cloud
)
:
    HeatTransferModel<CloudType>(dict, cloud, typeName),
    a_(this->coeffDict().template getOrDefault<scalar>("a", 2.0)),
    b_(this->coeffDict().template getOrDefault<scalar>("b", 0.6)),
    m_(this->coeffDict().template getOrDefault<scalar>("m", 1.0/2.0)),
    n_(this->coeffDict().template getOrDefault<scalar>("n", 1.0/3.0))
{}


template<class CloudType>
Foam::RanzMarshall<CloudType>::RanzMarshall(const RanzMarshall<CloudType>& htm)
:
    HeatTransferModel<CloudType>(htm),
    a_(htm.a_),
    b_(htm.b_),
    m_(htm.m_),
    n_(htm.n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::RanzMarshall<CloudType>::Nu
(
    const scalar Re,
    const scalar Pr
) const
{
    // (AOB:p. 18 below Eq. 42)
    return a_ + b_*pow(Re, m_)*pow(Pr, n_);
}


// ************************************************************************* //
