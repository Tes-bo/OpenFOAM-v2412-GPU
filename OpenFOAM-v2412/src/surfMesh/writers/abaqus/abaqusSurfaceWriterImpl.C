/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "IOmanip.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::surfaceWriters::abaqusWriter::writeFaceValue
(
    Ostream& os,
    const Type& value,
    const label elemId
) const
{
    if (fileFormats::ABAQUSCore::isEncodedSolidId(elemId))
    {
        const label solidId =
            fileFormats::ABAQUSCore::decodeSolidElementId(elemId);

        const label sideNum =
            fileFormats::ABAQUSCore::decodeSolidSideNum(elemId);

        // Element-id: 0-based to 1-based
        // Side-id: always 1-based
        os  << (solidId + 1) << ", P" << sideNum;
    }
    else
    {
        // Element-id: 0-based to 1-based
        os  << (elemId + 1) << ", P";
    }

    os << ", ";
    if (pTraits<Type>::nComponents == 1)
    {
        os << value;
    }
    else
    {
        os << mag(value);
    }
    os  << nl;

    return os;
}


template<class Type>
Foam::fileName Foam::surfaceWriters::abaqusWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Field:
    // 1) rootdir/<TIME>/<field>_surfaceName.inp
    // 2) rootdir/<field>/surfaceName_<TIME>.inp

    fileName outputFile;

    switch (outputLayout_)
    {
        case outputLayoutType::BY_TIME:
        {
            outputFile = outputPath_;
            if (useTimeDir() && !timeName().empty())
            {
                // Splice in time-directory
                outputFile =
                    outputPath_.path() / timeName() / outputPath_.name();
            }

            // Append <field>_surfaceName
            outputFile /= fieldName + '_' + outputPath_.name();
            break;
        }
        case outputLayoutType::BY_FIELD:
        {
            outputFile = outputPath_ / fieldName / outputPath_.name();
            if (!timeName().empty())
            {
                // Append time information to file name
                outputFile += '_' + timeName();
            }
            break;
        }
    }
    outputFile.ext("inp");


    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< " to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        const auto& values = tfield();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        // const scalar timeValue(0);


        // Additional bookkeeping for decomposing non tri/quad
        labelList decompOffsets;
        DynamicList<face> decompFaces;


        OFstream os(outputFile);

        if (noGeometry_ || wroteGeom_)
        {
            // Geometry already written (or suppressed)
            // - still need decomposition information

            fileFormats::ABAQUSCore::faceDecomposition
            (
                surf.points(),
                surf.faces(),
                decompOffsets,
                decompFaces
            );
        }
        else
        {
            // Write geometry (separate file)

            OFstream osGeom(outputFile.lessExt().ext("abq"));
            writeGeometry(osGeom, surf, decompOffsets, decompFaces);
        }

        // Write field
        // *DLOAD
        // Element-based: // elemId, P, 1000
        // Surface-based: // elemId, P4, 1000

        os  << "**" << nl
            << "** field = " << fieldName << nl
            << "** type = " << pTraits<Type>::typeName << nl;

        if (useTimeDir() && !timeName().empty())
        {
            os  << "** time = " << timeName() << nl;
        }

        os  << "**" << nl
            << "*DLOAD" << nl;


        // Regular (undecomposed) faces
        const faceList&    faces = surf.faces();
        const labelList& elemIds = surf.faceIds();

        // Possible to use faceIds?
        const bool useOrigFaceIds =
        (
            elemIds.size() == faces.size()
         && decompFaces.empty()
        );


        label elemId = 0;

        if (this->isPointData())
        {
            forAll(faces, facei)
            {
                if (useOrigFaceIds)
                {
                    // When available and not decomposed
                    elemId = elemIds[facei];
                }

                const label beginElemId = elemId;

                // Any face decomposition
                for
                (
                    label decompi = decompOffsets[facei];
                    decompi < decompOffsets[facei+1];
                    ++decompi
                )
                {
                    const face& f = decompFaces[decompi];

                    Type v = Zero;
                    for (const label verti : f)
                    {
                        v += values[verti];
                    }
                    v /= f.size();

                    writeFaceValue(os, v, elemId);  // 0-based
                    ++elemId;
                }


                // Face not decomposed
                if (beginElemId == elemId)
                {
                    const face& f = faces[facei];

                    Type v = Zero;
                    for (const label verti : f)
                    {
                        v += values[verti];
                    }
                    v /= f.size();

                    writeFaceValue(os, v, elemId);  // 0-based
                    ++elemId;
                }
            }
        }
        else
        {
            auto valIter = values.cbegin();

            forAll(faces, facei)
            {
                if (useOrigFaceIds)
                {
                    // When available and not decomposed
                    elemId = elemIds[facei];
                }

                const Type v(*valIter);
                ++valIter;

                label nValues =
                    max
                    (
                        label(1),
                        (decompOffsets[facei+1] - decompOffsets[facei])
                    );

                while (nValues--)
                {
                    writeFaceValue(os, v, elemId);  // 0-based
                    ++elemId;
                }
            }
        }

        os  << "**" << nl
            << "**" << nl;
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //
