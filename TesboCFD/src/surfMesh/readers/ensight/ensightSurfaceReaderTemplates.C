/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "SpanStream.H"
#include "ensightPTraits.H"

// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

template<class Type>
void Foam::ensightSurfaceReader::readFrom
(
    const std::string& buffer,
    Type& value
)
{
    ISpanStream is(buffer.data(), buffer.size());
    is  >> value;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::ensightSurfaceReader::readField
(
    const fileName& dataFile,
    const word& fieldName,
    const label timeIndex
) const
{
    auto tfield = tmp<Field<Type>>::New(surfPtr_->nFaces(), Zero);
    auto& field = tfield.ref();

    if (!masterOnly_ || UPstream::master(UPstream::worldComm))
    {
        // Use previously detected ascii/binary format
        ensightReadFile is(dataFile, readFormat_);

        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << is.name()
                << " for field " << fieldName
                << exit(FatalError);
        }

        // If transient single-file
        is.seekTime(timeIndex);


        // Check that data type is as expected
        // (assuming OpenFOAM generated the data set)
        string primitiveType;
        is.read(primitiveType);

        DebugInfo << "primitiveType: " << primitiveType << endl;

        if
        (
            debug
         && primitiveType != ensightPTraits<Type>::typeName
         && primitiveType != pTraits<Type>::typeName
        )
        {
            IOWarningInFunction(is)
                << "Expected <" << ensightPTraits<Type>::typeName
                << "> values for <" << pTraits<Type>::typeName
                << "> but found " << primitiveType << nl
                << "    This may be okay, but could indicate an error"
                << nl << nl;
        }

        string strValue;
        label intValue;

        // Read header info: part index, e.g. part 1
        is.read(strValue);
        is.read(intValue);

        label begFace = 0;

        // Loop through different element types when reading the field values
        for (const faceInfoTuple& facesInfo : faceTypeInfo_)
        {
            // [faceType, faceCount]
            const label endFace = begFace + facesInfo.second();

            DebugInfo
                << "Reading <" << pTraits<Type>::typeName << "> face type "
                << ensightFaces::elemNames[facesInfo.first()]
                << " data:" << facesInfo.second() << endl;

            if (begFace < endFace)
            {
                // The element type, optionally with 'undef'
                is.read(strValue);

                if (strValue.contains("undef"))
                {
                    // Skip undef entry
                    scalar value;
                    is.read(value);
                }

                // Ensight fields are written component-wise
                // (can be in different order than OpenFOAM uses)

                for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
                {
                    const direction cmpt =
                        ensightPTraits<Type>::componentOrder[d];

                    for (label facei = begFace; facei < endFace; ++facei)
                    {
                        scalar value;
                        is.read(value);
                        setComponent(field[facei], cmpt) = value;
                    }
                }

                begFace = endFace;
            }
        }
    }

    if (masterOnly_ && UPstream::parRun())
    {
        Pstream::broadcast(field, UPstream::worldComm);
    }

    return tfield;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::ensightSurfaceReader::readField
(
    const label timeIndex,
    const label fieldIndex
) const
{
    if (fieldIndex < 0 || fieldIndex >= fieldNames_.size())
    {
        FatalErrorInFunction
            << "Invalid timeIndex:" << timeIndex
            << " should be in range [0.." << fieldNames_.size() << ')' << nl
            << "Possibly used incorrect field lookup name. Known field names: "
            << flatOutput(fieldNames_) << nl
            << exit(FatalError);
    }

    const word& fieldName = fieldNames_[fieldIndex];

    const label fileIndex =
    (
        (timeIndex >= 0 && timeIndex < fileNumbers_.size())
      ? fileNumbers_[timeIndex]
      : (timeStartIndex_ + timeIndex*timeIncrement_)
    );

    const fileName dataFile
    (
        baseDir_
      / ensightCase::expand_mask(fieldFileNames_[fieldIndex], fileIndex)
    );

    if (debug)
    {
        Pout<< "Read <" << pTraits<Type>::typeName << "> field, file="
            << dataFile << endl;
    }

    return readField<Type>(dataFile, fieldName, timeIndex);
}


// ************************************************************************* //
