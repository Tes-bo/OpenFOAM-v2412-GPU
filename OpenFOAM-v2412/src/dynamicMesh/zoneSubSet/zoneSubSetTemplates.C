/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "volFieldsFwd.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::Detail::zoneSubSet::mapToZone
(
    const GeometricField<Type, fvPatchField, volMesh>& subVolField
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    // The full-mesh field, with zero for non-zoned areas
    auto tcellZonesField = volFieldType::New
    (
        subVolField.name(),
        subsetter_.baseMesh(),
        dimensioned<Type>(subVolField.dimensions())
    );
    auto& cellZonesField = tcellZonesField.ref();


    // Map field in global mesh on original zones
    // - without any halo cells

    const labelList& cellMap = subsetter_.cellMap();

    forAll(cellMap, subCelli)
    {
        const label celli = cellMap[subCelli];

        if (!haloCells_.test(celli))
        {
            cellZonesField[celli] = subVolField[subCelli];
        }
    }

    return tcellZonesField;
}


// ************************************************************************* //
