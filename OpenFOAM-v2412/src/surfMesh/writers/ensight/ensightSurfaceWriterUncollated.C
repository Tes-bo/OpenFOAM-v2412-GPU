/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::ensightWriter::writeUncollated()
{
    checkOpen();

    const ensight::FileName baseName(outputPath_.name());


    // Uncollated
    // ==========
    // CaseFile:  rootdir/<TIME>/NAME.case
    // Geometry:  rootdir/<TIME>/NAME.00000000.mesh

    fileName outputDir;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputDir = outputPath_.path() / timeName();
    }
    else
    {
        outputDir = outputPath_.path();
    }

    const fileName outputFile = outputDir / baseName + ".case";

    if (verbose_)
    {
        Info<< "Writing case file to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!Foam::isDir(outputDir))
        {
            Foam::mkDir(outputDir);
        }

        // The geometry
        ensightOutputSurface part
        (
            surf.points(),
            surf.faces(),
            baseName
        );

        // Two-argument form for path-name to avoid validating outputDir
        ensightGeoFile osGeom
        (
            outputDir,
            baseName + ".00000000.mesh",
            caseOpts_.format()
        );

        osGeom.beginGeometry();
        part.write(osGeom);  // serial

        // Update case file
        OFstream osCase
        (
            IOstreamOption::ATOMIC,
            outputFile,
            IOstreamOption::ASCII
        );
        ensightCase::setTimeFormat(osCase, caseOpts_);  // time-format

        osCase
            << "FORMAT" << nl
            << "type: ensight gold" << nl
            << nl
            << "GEOMETRY" << nl
            << "model:        1     " << osGeom.name().name() << nl
            << nl
            << "TIME" << nl;

        ensightCase::printTimeset(osCase, 1, scalar(0));
    }

    wroteGeom_ = true;
    return outputFile;
}


template<class Type>
Foam::fileName Foam::surfaceWriters::ensightWriter::writeUncollated
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    const ensight::FileName baseName(outputPath_.name());
    const ensight::VarName  varName(fieldName);


    // Uncollated
    // ==========
    // CaseFile:  rootdir/time/<field>/NAME.case
    // Geometry:  rootdir/time/<field>/NAME.<index>.mesh
    // Field:     rootdir/time/<field>/NAME.<index>.<field>

    // Variable name as sub-directory for results. Eg,
    // - VAR1/NAME1.case
    // - VAR1/NAME1.00000000.mesh
    // - VAR1/NAME1.00000001.VAR1
    // and
    // - VAR2/NAME1.case
    // - VAR2/NAME1.00000000.mesh
    // - VAR2/NAME1.00000001.VAR2

    fileName outputDir;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputDir = outputPath_.path() / timeName();
    }
    else
    {
        outputDir = outputPath_.path();
    }

    const fileName baseDir = outputDir / varName;
    const word   timeDir = timeName();
    const scalar timeValue = currTime_.value();

    const fileName outputFile = baseDir / baseName + ".case";

    if (verbose_)
    {
        Info<< "Writing case file to " << outputFile << nl;
    }

    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< endl;
    }

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!Foam::isDir(outputFile.path()))
        {
            Foam::mkDir(outputFile.path());
        }

        // The geometry
        ensightOutputSurface part
        (
            surf.points(),
            surf.faces(),
            baseName
        );

        // Two-argument form for path-name to avoid validating base-dir
        ensightGeoFile osGeom
        (
            baseDir,
            baseName + ".00000000.mesh",
            caseOpts_.format()
        );
        ensightFile osField
        (
            baseDir,
            baseName + ".00000000." + varName,
            caseOpts_.format()
        );

        osGeom.beginGeometry();
        part.write(osGeom);  // serial

        // Write field (serial)
        osField.write(ensightPTraits<Type>::typeName);
        osField.newline();
        part.writeData(osField, tfield(), this->isPointData());


        // Update case file
        {
            OFstream osCase
            (
                IOstreamOption::ATOMIC,
                outputFile,
                IOstreamOption::ASCII
            );
            ensightCase::setTimeFormat(osCase, caseOpts_);  // time-format

            osCase
                << "FORMAT" << nl
                << "type: ensight gold" << nl
                << nl
                << "GEOMETRY" << nl
                << "model:  1   " << osGeom.name().name() << nl
                << nl
                << "VARIABLE" << nl
                << ensightPTraits<Type>::typeName
                <<
                (
                    this->isPointData()
                  ? " per node:    1  "  // time-set 1
                  : " per element: 1  "  // time-set 1
                )
                << setw(15) << varName << ' '
                << baseName.c_str() << ".********."
                << ensight::FileName(varName).c_str() << nl;

            osCase
                << nl
                << "TIME" << nl;

                ensightCase::printTimeset(osCase, 1, timeValue);
                osCase << "# end" << nl;
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //
