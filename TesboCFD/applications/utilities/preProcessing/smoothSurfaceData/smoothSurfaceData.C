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

Application
    smoothSurfaceData

Description
    Pre-processing, filtering of surface field data.
    Currently this is just used to test filter settings
    (for MappedFile) and only generates vtk output, which can be
    easily loaded in ParaView.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "clockTime.H"
#include "primitiveFields.H"
#include "surfaceReader.H"
#include "surfaceWriter.H"
#include "foamVtkSurfaceWriter.H"
#include "MappedFileFilterField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Testing, pre-processing, filtering of surface field data"
    );

    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();

    argList::addArgument("input", "The input surface file");

    argList::addOption("radius", "m", "Specify filter radius [metres]");
    argList::addOption("sweeps", "N", "Number of median filter stages");
    argList::addOption
    (
        "field",
        "name",
        "Field <scalar> to process (default: T)"
    );

    argList::addOption
    (
        "read-format",
        "type",
        "Input format (default: ensight)"
    );

    argList args(argc, argv);
    // #include "setRootCase.H"

    const int optVerbose = args.verbose();

    const word readFileType
    (
        args.getOrDefault<word>("read-format", "ensight")
    );

    // Constant radius searching
    label filterSweeps_(1);
    scalar filterRadius_(0);
    args.readIfPresent("sweeps", filterSweeps_);
    args.readIfPresent("radius", filterRadius_);

    Info<< nl
        << "Filter: radius=" << filterRadius_
        << " sweeps=" << filterSweeps_ << endl;

    // Simple sanity check
    if ((filterSweeps_ < 1) || (filterRadius_ <= VSMALL))
    {
        Info<< nl << "Nothing to do. Exiting..." << nl << endl;
        return 0;
    }

    word fieldName("T");
    args.readIfPresent("field", fieldName);


    const fileName importName = args.get<fileName>(1);

    auto readerPtr_ = surfaceReader::New(readFileType, importName);

    auto& reader = readerPtr_();

    const label fieldIndex = reader.fieldNames(0).find(fieldName);
    if (fieldIndex == -1)
    {
        FatalErrorInFunction
            << "Unable to find field name: " << fieldName
            << " in list of available fields: " << reader.fieldNames(0)
            << exit(FatalError);
    }

    clockTime timing;

    const meshedSurface& geom = reader.geometry(0);

    Info<< nl << "Read " << geom.nFaces() << " faces and "
        << geom.nPoints() << " points in "
        << timing.timeIncrement() << "s" << endl;

    PatchFunction1Types::FilterField fieldFilter;

    PatchFunction1Types::FilterField::debug = optVerbose;

    fieldFilter.reset(geom, filterRadius_);

    Info<< nl << "Built weights/addressing "
        << timing.timeIncrement() << "s" << endl;

    Info<< nl << "Processing " << reader.times().size() << " times" << nl;

    const instantList times(reader.times());

    forAll(times, timeIndex)
    {
        tmp<scalarField> tfield
        (
            reader.field(timeIndex, fieldIndex, pTraits<scalar>::zero)
        );

        tfield = fieldFilter.evaluate(tfield, filterSweeps_);

        vtk::surfaceWriter writer
        (
            geom.points(),
            geom,
            word::printf("filtered_%06d", timeIndex)
        );

        writer.piece(geom.points(), geom);

        writer.writeGeometry();
        writer.beginCellData();
        writer.writeCellData(fieldName, tfield());
    }

    Info<< nl << "Smoothing/writing "
            << timing.timeIncrement() << "s" << endl;

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
