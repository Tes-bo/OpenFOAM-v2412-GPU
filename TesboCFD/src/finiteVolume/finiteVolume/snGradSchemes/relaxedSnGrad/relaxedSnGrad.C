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

#include "fv.H"
#include "relaxedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::relaxedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfFieldType;

    // Calculate explicit correction field
    tmp<SurfFieldType> tcorrection = correctedScheme_().correction(vf);

    // Retrieve relaxation factor value
    const word fieldName(vf.name());
    const word oldFieldName(fieldName + "_0");
    const scalar relax =
        vf.mesh().fieldRelaxationFactor("snGrad("+fieldName+")");

    // Return explicit correction field if
    // previous-time step correction is unavailable
    const objectRegistry& obr = vf.db();
    if (!obr.foundObject<SurfFieldType>(oldFieldName))
    {
        SurfFieldType* oldCorrection =
            new SurfFieldType(oldFieldName, tcorrection());
        oldCorrection->store();
    }

    // Return under/over-relaxed explicit correction field
    auto trelaxedCorrection = tmp<SurfFieldType>::New(tcorrection());

    SurfFieldType& oldCorrection =
        obr.lookupObjectRef<SurfFieldType>(oldFieldName);

    trelaxedCorrection.ref() *= relax;
    trelaxedCorrection.ref() += (scalar(1) - relax)*oldCorrection;

    oldCorrection = tcorrection;

    return trelaxedCorrection;
}


// ************************************************************************* //
