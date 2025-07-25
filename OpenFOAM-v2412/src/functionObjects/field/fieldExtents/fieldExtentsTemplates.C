/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "fieldExtents.H"
#include "volFields.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField> Foam::functionObjects::fieldExtents::calcMask
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return
        pos
        (
            mag(field)
          - dimensionedScalar("t", field.dimensions(), threshold_)
        );
}


template<class Type>
void Foam::functionObjects::fieldExtents::calcFieldExtents
(
    const word& fieldName,
    const bool calcMag
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* fieldPtr = obr_.findObject<VolFieldType>(fieldName);

    if (!fieldPtr)
    {
        return;
    }

    auto extents = [this](const scalarField& mask, const vectorField& C)
    {
        boundBox extents;
        forAll(mask, i)
        {
            if (mask[i] > 0.5)
            {
                extents.add(C[i] - C0_);
            }
        };

        extents.reduce();

        if (extents.empty())
        {
            extents.reset(point::zero);
        }

        return extents;
    };

    Log << "field: " << fieldName << nl;

    writeCurrentTime(file());

    tmp<volScalarField> tmask = calcMask<Type>(*fieldPtr);
    const volScalarField& mask = tmask();

    // Internal field
    if (internalField_)
    {
        boundBox bb(extents(mask, mesh_.C()));
        Log << "    internal field: " << bb << nl;
        file() << bb;

        this->setResult(fieldName + "_internal_min", bb.min());
        this->setResult(fieldName + "_internal_max", bb.max());
    }

    // Patches
    for (const label patchi : patchIDs_)
    {
        const fvPatchScalarField& maskp = mask.boundaryField()[patchi];
        boundBox bb(extents(maskp, maskp.patch().Cf()));
        const word& patchName = maskp.patch().name();
        Log << "    patch " << patchName << ": " << bb << nl;
        file() << bb;
        this->setResult(fieldName + "_" + patchName + "_min", bb.min());
        this->setResult(fieldName + "_" + patchName + "_max", bb.max());
    }

    Log << endl;
    file() << endl;
}


// ************************************************************************* //
