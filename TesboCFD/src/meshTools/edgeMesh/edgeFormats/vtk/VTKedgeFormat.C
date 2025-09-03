/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "VTKedgeFormat.H"
#include "Fstream.H"
#include "Time.H"
#include "foamVtkLineWriter.H"
#include "vtkUnstructuredReader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::VTKedgeFormat::VTKedgeFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::VTKedgeFormat::read
(
    const fileName& filename
)
{
    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // Use dummy Time for objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New());

    objectRegistry obr
    (
        IOobject
        (
            "vtk::edgeFormat",
            *dummyTimePtr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // Construct reader to read file
    vtkUnstructuredReader reader(obr, is);


    // Extract lines
    storedPoints().transfer(reader.points());

    label nEdges = 0;
    for (const auto& lineVerts : reader.lines())
    {
        if (lineVerts.size() > 1)
        {
            nEdges += (lineVerts.size()-1);
        }
    }
    storedEdges().resize(nEdges);

    nEdges = 0;
    for (const auto& lineVerts : reader.lines())
    {
        for (label i = 1; i < lineVerts.size(); ++i)
        {
            storedEdges()[nEdges++] = edge(lineVerts[i-1], lineVerts[i]);
        }
    }

    return true;
}


void Foam::fileFormats::VTKedgeFormat::write
(
    const fileName& filename,
    const edgeMesh& eMesh,
    IOstreamOption,
    const dictionary& options
)
{
    // NB: restrict output to legacy ascii so that we are still able
    // to read it with vtkUnstructuredReader

    vtk::outputOptions opts(vtk::formatType::LEGACY_ASCII);

    opts.precision
    (
        options.getOrDefault("precision", IOstream::defaultPrecision())
    );

    vtk::lineWriter writer
    (
        eMesh.points(),
        eMesh.edges(),
        opts,
        filename,
        false  // non-parallel write (edgeMesh already serialized)
    );

    writer.beginFile("OpenFOAM edgeMesh");
    writer.writeGeometry();
}


// ************************************************************************* //
