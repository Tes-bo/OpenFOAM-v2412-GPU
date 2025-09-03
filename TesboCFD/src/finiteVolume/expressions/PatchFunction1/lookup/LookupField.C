/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "volFields.H"
#include "volMesh.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::LookupField<Type>::LookupField
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    lookupBase(dict)
{}


template<class Type>
Foam::PatchFunction1Types::LookupField<Type>::LookupField
(
    const LookupField<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    lookupBase(rhs)
{}


template<class Type>
Foam::PatchFunction1Types::LookupField<Type>::LookupField
(
    const LookupField<Type>& rhs
)
:
    LookupField<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::LookupField<Type>::value(const scalar x) const
{
    typedef UniformDimensionedField<Type> UType;

    const objectRegistry& db = patchFunction1Base::obr();
    const label patchi = patchFunction1Base::patch().index();

    if (this->faceValues())
    {
        typedef GeometricField<Type, fvPatchField, volMesh> VType;
        typedef GeometricField<Type, fvsPatchField, surfaceMesh> SType;

        // Try:
        // - as volField in local scope
        // - as surfaceField in local scope
        // - as UniformDimensionedField recursively

        const regIOobject* ptr = db.cfindIOobject(name_, false);

        if (ptr)
        {
            const auto* vPtr = dynamic_cast<const VType*>(ptr);
            if (vPtr)
            {
                return vPtr->boundaryField()[patchi];
            }

            const auto* sPtr = dynamic_cast<const SType*>(ptr);
            if (sPtr)
            {
                return sPtr->boundaryField()[patchi];
            }

            const auto* uPtr = dynamic_cast<const UType*>(ptr);
            if (uPtr)
            {
                return Field<Type>(this->size(), uPtr->value());
            }
        }

        // Done db level. Try recursion
        ptr = db.parent().cfindIOobject(name_, true);

        if (ptr)
        {
            const auto* uPtr = dynamic_cast<const UType*>(ptr);
            if (uPtr)
            {
                return Field<Type>(this->size(), uPtr->value());
            }
        }

        FatalErrorInFunction
            << nl
            << "    failed lookup of " << name_
            << " (objectRegistry "
            << db.name()
            << ")\n    available objects of type " << VType::typeName
            << ':' << nl
            << db.names<VType>() << nl
            << "    available objects of type " << SType::typeName
            << ':' << nl
            << db.names<SType>() << nl
            << "    available objects of type " << UType::typeName
            << ':' << nl
            << db.names<UType>() << nl
            << exit(FatalError);
        return Field<Type>::null();
    }
    else
    {
        // Assume pointField
        typedef GeometricField<Type, pointPatchField, pointMesh> PType;

        const regIOobject* ptr = db.cfindIOobject(name_, false);

        if (ptr)
        {
            const auto* pPtr = dynamic_cast<const PType*>(ptr);
            if (pPtr)
            {
                return pPtr->boundaryField()[patchi].patchInternalField();
            }
        }

        // Re-do as uniform field. Note: could repeat logic above
        const auto& obj = db.lookupObject<UType>(name_, true);

        return Field<Type>(this->size(), obj.value());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::LookupField<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return (x2-x1)*value(0.5*(x1+x2));
}


template<class Type>
void Foam::PatchFunction1Types::LookupField<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);
    os.writeEntry(this->name(), type());
    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
