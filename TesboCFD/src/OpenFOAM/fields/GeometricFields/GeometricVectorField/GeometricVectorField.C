/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "GeometricVectorField.H"
#include "vectorFieldField.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::zip
(
    GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& result,
    const GeometricField<Cmpt, PatchField, GeoMesh>& x,
    const GeometricField<Cmpt, PatchField, GeoMesh>& y,
    const GeometricField<Cmpt, PatchField, GeoMesh>& z
)
{
    Foam::zip
    (
        result.primitiveFieldRef(),
        x.primitiveField(),
        y.primitiveField(),
        z.primitiveField()
    );

    Foam::zip
    (
        result.boundaryFieldRef(),
        x.boundaryField(),
        y.boundaryField(),
        z.boundaryField()
    );
}


template<class Cmpt, template<class> class PatchField, class GeoMesh>
void Foam::unzip
(
    const GeometricField<Vector<Cmpt>, PatchField, GeoMesh>& input,
    GeometricField<Cmpt, PatchField, GeoMesh>& x,
    GeometricField<Cmpt, PatchField, GeoMesh>& y,
    GeometricField<Cmpt, PatchField, GeoMesh>& z
)
{
    Foam::unzip
    (
        input.primitiveField(),
        x.primitiveFieldRef(),
        y.primitiveFieldRef(),
        z.primitiveFieldRef()
    );

    Foam::unzip
    (
        input.boundaryField(),
        x.boundaryFieldRef(),
        y.boundaryFieldRef(),
        z.boundaryFieldRef()
    );
}


// ************************************************************************* //
