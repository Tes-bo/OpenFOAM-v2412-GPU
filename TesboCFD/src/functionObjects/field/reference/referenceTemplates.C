/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "Function1.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::reference::calcType()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* vfPtr = findObject<VolFieldType>(fieldName_);

    if (vfPtr)
    {
        const VolFieldType& vf = *vfPtr;

        // non-mandatory
        dimensioned<Type> offset("offset", vf.dimensions(), Zero, localDict_);

        dimensioned<Type> refValue("refValue", vf.dimensions(), Zero);

        autoPtr<Function1<Type>> valuePtr
        (
            Function1<Type>::New("refValue", localDict_, &this->mesh_)
        );

        refValue.value() = valuePtr->value(this->time().value());

        Info<< "    Reference value: " << refValue.value() << endl;

        return store
        (
            resultName_,
            scale_*(vf - refValue + offset)
        );
    }

    return false;
}


// ************************************************************************* //
