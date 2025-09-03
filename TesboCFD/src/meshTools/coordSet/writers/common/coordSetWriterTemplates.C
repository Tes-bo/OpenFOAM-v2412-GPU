/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "transformField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coordSetWriter::adjustFieldTemplate
(
    const word& fieldName,
    const tmp<Field<Type>>& tfield
) const
{
    if (verbose_)
    {
        Info<< "Writing field " << fieldName;
    }

    tmp<Field<Type>> tadjusted;

    // Output scaling for the variable, but not for integer types
    // which are typically ids etc.
    if (!std::is_integral<Type>::value)
    {
        scalar value;

        // Remove *uniform* reference level
        if
        (
            fieldLevel_.readIfPresent(fieldName, value, keyType::REGEX)
         && !equal(value, 0)
        )
        {
            // Could also detect brackets (...) and read accordingly
            // or automatically scale by 1/sqrt(nComponents) instead ...

            Type refLevel;
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
            {
                setComponent(refLevel, cmpt) = value;
            }

            if (verbose_)
            {
                Info<< " [level " << refLevel << ']';
            }

            if (!tadjusted)
            {
                // Steal or clone
                tadjusted.reset(tfield.ptr());
            }

            // Remove offset level
            tadjusted.ref() -= refLevel;
        }

        // Apply scaling
        if
        (
            fieldScale_.readIfPresent(fieldName, value, keyType::REGEX)
         && !equal(value, 1)
        )
        {
            if (verbose_)
            {
                Info<< " [scaling " << value << ']';
            }

            if (!tadjusted)
            {
                // Steal or clone
                tadjusted.reset(tfield.ptr());
            }

            // Apply scaling
            tadjusted.ref() *= value;
        }

        // Rotate fields (vector and non-spherical tensors)
        if
        (
            (is_vectorspace<Type>::value && pTraits<Type>::nComponents > 1)
         && geometryTransform_.good()
         && !geometryTransform_.R().is_identity()
        )
        {
            if (!tadjusted)
            {
                // Steal or clone
                tadjusted.reset(tfield.ptr());
            }

           Foam::transform
           (
               tadjusted.ref(),
               geometryTransform_.R(),
               tadjusted()
           );
        }
    }

    return (tadjusted ? tadjusted : tfield);
}


template<class Type>
Foam::UPtrList<const Foam::Field<Type>>
Foam::coordSetWriter::repackageFields(const Field<Type>& field)
{
    UPtrList<const Field<Type>> fieldPtrs(1);
    fieldPtrs.set(0, &field);

    return fieldPtrs;
}


template<class Type>
Foam::UPtrList<const Foam::Field<Type>>
Foam::coordSetWriter::repackageFields(const UList<Field<Type>>& fieldValues)
{
    UPtrList<const Field<Type>> fieldPtrs(fieldValues.size());
    forAll(fieldValues, i)
    {
        fieldPtrs.set(i, &(fieldValues[i]));
    }

    return fieldPtrs;
}


template<class Type>
void Foam::coordSetWriter::writeTable
(
    Ostream& os,
    const coordSet& coords,
    const UList<Type>& values,
    const char* sep
)
{
    forAll(coords, pointi)
    {
        // Output coordinate (point or scalar) with separator
        if (coords.hasVectorAxis())
        {
            const vector& p = coords.vectorCoord(pointi);
            os << p.x() << sep << p.y() << sep << p.z();
        }
        else
        {
            os << coords.scalarCoord(pointi);
        }

        // Output component values with separator
        const auto& val = values[pointi];
        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            os << sep << component(val, d);
        }
        os << nl;
    }
}


// ************************************************************************* //
