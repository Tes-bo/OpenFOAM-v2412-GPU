/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

inline Foam::label Foam::expressions::patchExpr::parseDriver::size
(
    const FieldAssociation geoType
) const
{
    switch (geoType)
    {
        case FieldAssociation::POINT_DATA :
            return patch_.patch().nPoints();

        case FieldAssociation::FACE_DATA :
            return patch_.patch().size();

        default:
            break;
    }
    return 0;
}


inline Foam::tmp<Foam::boolField>
Foam::expressions::patchExpr::parseDriver::parseDriver::field_cellSet
(
    const word& name
) const
{
    return field_cellSelection
    (
        name,
        topoSetSource::sourceType::CELLSET_SOURCE
    );
}


inline Foam::tmp<Foam::boolField>
Foam::expressions::patchExpr::parseDriver::field_cellZone
(
    const word& name
) const
{
    return field_cellSelection
    (
        name,
        topoSetSource::sourceType::CELLZONE_SOURCE
    );
}


inline Foam::tmp<Foam::boolField>
Foam::expressions::patchExpr::parseDriver::field_faceSet
(
    const word& name
) const
{
    return field_faceSelection
    (
        name,
        topoSetSource::sourceType::FACESET_SOURCE
    );
}


inline Foam::tmp<Foam::boolField>
Foam::expressions::patchExpr::parseDriver::field_faceZone
(
    const word& name
) const
{
    return field_faceSelection
    (
        name,
        topoSetSource::sourceType::FACEZONE_SOURCE
    );
}


// ************************************************************************* //
