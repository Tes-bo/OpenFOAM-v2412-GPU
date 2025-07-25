/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "indexedVertexEnum.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::indexedVertexEnum::vertexType
>
Foam::indexedVertexEnum::vertexTypeNames_
({
    { vertexType::vtUnassigned, "Unassigned" },
    { vertexType::vtInternal, "Internal" },
    { vertexType::vtInternalNearBoundary, "InternalNearBoundary" },
    { vertexType::vtInternalSurface, "InternalSurface" },
    { vertexType::vtInternalSurfaceBaffle, "InternalSurfaceBaffle" },
    { vertexType::vtExternalSurfaceBaffle, "ExternalSurfaceBaffle" },
    { vertexType::vtInternalFeatureEdge, "InternalFeatureEdge" },
    { vertexType::vtInternalFeatureEdgeBaffle, "InternalFeatureEdgeBaffle" },
    { vertexType::vtExternalFeatureEdgeBaffle, "ExternalFeatureEdgeBaffle" },
    { vertexType::vtInternalFeaturePoint, "InternalFeaturePoint" },
    { vertexType::vtExternalSurface, "ExternalSurface" },
    { vertexType::vtExternalFeatureEdge, "ExternalFeatureEdge" },
    { vertexType::vtExternalFeaturePoint, "ExternalFeaturePoint" },
    { vertexType::vtFar, "Far" },
    { vertexType::vtConstrained, "Constrained" },
});


const Foam::Enum
<
    Foam::indexedVertexEnum::vertexMotion
>
Foam::indexedVertexEnum::vertexMotionNames_
{
    { vertexMotion::fixed, "fixed" },
    { vertexMotion::movable, "movable" },
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Foam::indexedVertexEnum::vertexType& v
)
{
    os  << static_cast<int>(v);

    return os;
}

Foam::Istream& Foam::operator>>
(
    Istream& is,
    Foam::indexedVertexEnum::vertexType& v
)
{
    int type;
    is  >> type;

    v = static_cast<Foam::indexedVertexEnum::vertexType>(type);

    return is;
}

// ************************************************************************* //
