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

#include "IOField.H"
#include "ensightOutputCloud.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::functionObjects::ensightCloudWriteObject::writeFields
(
    const word& cloudName,
    const objectRegistry& obrTmp
)
{
    static_assert
    (
        (
            std::is_same<label, typename pTraits<Type>::cmptType>::value
         || std::is_floating_point<typename pTraits<Type>::cmptType>::value
        ),
        "Label and Floating-point vector space only"
    );

    // Other integral types (eg, bool etc) would need cast/convert to label.
    // Similarly for labelVector etc.


    // Fields are not always on all processors (eg, multi-component parcels).
    // Thus need to resolve names between all processors.

    wordList fieldNames =
    (
        selectFields_.size()
      ? obrTmp.names<IOField<Type>>(selectFields_)
      : obrTmp.names<IOField<Type>>()
    );

    Pstream::combineReduce(fieldNames, ListOps::uniqueEqOp<word>());
    Foam::sort(fieldNames);  // Consistent order

    DynamicList<Type> scratch;

    for (const word& fieldName : fieldNames)
    {
        const List<Type>* fldPtr = obrTmp.findObject<IOField<Type>>(fieldName);
        const List<Type>& values = (fldPtr ? *fldPtr : List<Type>::null());

        autoPtr<ensightFile> os =
            ensCase().newCloudData<Type>(cloudName, fieldName);

        if (applyFilter_)
        {
            scratch.resize_nocopy(parcelAddr_.count());

            auto iter = scratch.begin();

            for (const label idx : parcelAddr_)
            {
                *iter = values[idx];
                ++iter;
            }

            // TBD:
            // recalculate globalIndex instead of relying on procAddr_ ?

            ensightOutput::writeCloudField(os.ref(), scratch, procAddr_);
        }
        else
        {
            // TBD:
            // recalculate globalIndex instead of relying on procAddr_ ?

            ensightOutput::writeCloudField(os.ref(), values, procAddr_);
        }
    }

    return fieldNames;
}


// ************************************************************************* //
