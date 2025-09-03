/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include <fstream>
#include "foamVtkWriteTopoSet.H"
#include "polyMesh.H"
#include "pointSet.H"
#include "globalIndex.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::vtk::writePointSet
(
    const polyMesh& mesh,
    const pointSet& set,
    const vtk::outputOptions opts,
    const fileName& file,
    bool parallel
)
{
    vtk::outputOptions opts_(opts);
    opts_.append(false);  // Do not support append

    const bool legacy = opts_.legacy();

    // Only allow parallel if really is a parallel run.
    parallel = parallel && Pstream::parRun();


    std::ofstream os_;
    autoPtr<vtk::formatter> format;

    // Open a file and attach a formatter
    // - on master (always)
    // - on slave if not parallel
    //
    // This means we can always check if format_ is defined to know if output
    // is desired on any particular process.

    if (Pstream::master() || !parallel)
    {
        mkDir(file.path());

        // Extension is inappropriate. Legacy instead of xml, or vice versa.
        const word ext = vtk::fileExtension[vtk::fileTag::POLY_DATA];

        if (file.has_ext(ext))
        {
            // Extension is correct
            os_.open(file);
        }
        else if
        (
            legacy
          ? file.has_ext(ext)
          : file.has_ext(vtk::legacy::fileExtension)
        )
        {
            // Extension is inappropriate. Legacy instead of xml, or vice versa.
            os_.open(file.lessExt() + "." + ext);
        }
        else
        {
            // Extension added automatically
            os_.open(file + "." + ext);
        }

        format = opts_.newFormatter(os_);
    }


    //-------------------------------------------------------------------------

    labelField pointLabels(set.sortedToc());

    label numberOfPoints = pointLabels.size();

    if (parallel)
    {
        reduce(numberOfPoints, sumOp<label>());
    }

    if (format)
    {
        const auto& title = set.name();

        if (legacy)
        {
            // beginFile:

            legacy::fileHeader<vtk::fileTag::POLY_DATA>(format(), title);

            // beginPoints:

            legacy::beginPoints(os_, numberOfPoints);
        }
        else
        {
            // XML (inline)

            // beginFile:

            format()
                .xmlHeader()
                .xmlComment(title)
                .beginVTKFile<vtk::fileTag::POLY_DATA>();

            // beginPiece:
            format()
                .tag
                (
                    vtk::fileTag::PIECE,
                    vtk::fileAttr::NUMBER_OF_POINTS, numberOfPoints
                );

            // beginPoints:
            const uint64_t payLoad = vtk::sizeofData<float,3>(numberOfPoints);

            format()
                .tag(vtk::fileTag::POINTS)
                .beginDataArray<float,3>(vtk::dataArrayAttr::POINTS);
            format().writeSize(payLoad);
        }
    }


    //-------------------------------------------------------------------------

    // pointLabels are the addressing for an indirect list

    if (parallel)
    {
        vtk::writeListParallel(format.ref(), mesh.points(), pointLabels);
    }
    else
    {
        vtk::writeList(format(), mesh.points(), pointLabels);
    }

    if (format)
    {
        format().flush();
        format().endDataArray();

        if (!legacy)
        {
            format()
                .endTag(vtk::fileTag::POINTS);
        }


        // beginPointData:
        if (legacy)
        {
            legacy::beginPointData(format(), numberOfPoints, 1); // 1 field
        }
        else
        {
            format().beginPointData();
        }
    }

    if (format)
    {
        // pointID

        if (legacy)
        {
            // 1 component
            legacy::intField<1>(format(), "pointID", numberOfPoints);
        }
        else
        {
            const uint64_t payLoad = vtk::sizeofData<label>(numberOfPoints);

            format().beginDataArray<label>("pointID");
            format().writeSize(payLoad);
        }
    }


    if (parallel)
    {
        const globalIndex pointIdOffset(mesh.nPoints());

        vtk::writeListParallel(format.ref(), pointLabels, pointIdOffset);
    }
    else
    {
        vtk::writeList(format(), pointLabels);
    }


    if (format)
    {
        format().flush();
        format().endDataArray();
        format().endPointData();
        format().endPiece();

        format().endTag(vtk::fileTag::POLY_DATA)
            .endVTKFile();
    }

    return true;
}


// ************************************************************************* //
