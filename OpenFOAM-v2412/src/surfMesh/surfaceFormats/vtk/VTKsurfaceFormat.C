/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "VTKsurfaceFormat.H"
#include "vtkUnstructuredReader.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "faceTraits.H"
#include <fstream>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::writePolys
(
    vtk::formatter& format,
    const UList<Face>& faces
)
{
    // connectivity count without additional storage (done internally)
    label nConnectivity = 0;
    for (const Face& f : faces)
    {
        nConnectivity += f.size();
    }

    vtk::legacy::beginPolys
    (
        format.os(),
        faces.size(),
        nConnectivity
    );


    // legacy: size + connectivity together
    // [nPts, id1, id2, ..., nPts, id1, id2, ...]

    for (const Face& f : faces)
    {
        format.write(f.size());  // The size prefix
        vtk::writeList(format, f);
    }

    format.flush();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::VTKsurfaceFormat<Face>::VTKsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::VTKsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    // Clear everything
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename << nl
            << exit(FatalError);
    }

    // Assume groups are not intermixed
    bool sorted = true;


    // Use dummy Time for objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New());

    objectRegistry obr
    (
        IOobject
        (
            "vtk::surfaceFormat",
            *dummyTimePtr,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        )
    );

    // Read all
    vtkUnstructuredReader reader(obr, is);
    const faceList& faces = reader.faces();

    // Assume all faces in zone0 unless a region field is present
    labelList zones(faces.size(), Zero);

    for (auto fieldName : { "region", "STLSolidLabeling" })
    {
        const labelIOField* lptr =
            reader.cellData().findObject<labelIOField>(fieldName);

        if (lptr)
        {
            label i = 0;
            for (const auto& region : *lptr)
            {
                zones[i++] = label(region);
            }
            break;
        }

        const scalarIOField* sptr =
            reader.cellData().findObject<scalarIOField>(fieldName);

        if (sptr)
        {
            label i = 0;
            for (const auto& region : *sptr)
            {
                zones[i++] = label(region);
            }
            break;
        }
    }


    // Create zone names
    const label nZones = max(zones)+1;
    wordList zoneNames(nZones);
    forAll(zoneNames, i)
    {
        zoneNames[i] = surfZone::defaultName(i);
    }


    // Check if it needs triangulation
    label nTri = 0;
    if (faceTraits<Face>::isTri())
    {
        for (const face& f : faces)
        {
            nTri += f.nTriangles();
        }
    }

    DynamicList<label> dynElemId; // unused

    if (nTri > faces.size())
    {
        // We are here if the target surface needs triangles and
        // the source surface has non-triangles

        DynamicList<Face> dynFaces(nTri);
        DynamicList<label> dynZones(nTri);

        forAll(faces, facei)
        {
            const face& f = faces[facei];
            for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
            {
                const label fp2 = f.fcIndex(fp1);

                dynFaces.append(Face{f[0], f[fp1], f[fp2]});
                dynZones.append(zones[facei]);
            }
        }
        zones.clear();

        // Count
        labelList zoneSizes(nZones, Zero);
        for (const label zonei : dynZones)
        {
            zoneSizes[zonei]++;
        }

        this->sortFacesAndStore(dynFaces, dynZones, dynElemId, sorted);

        // Add zones (retaining empty ones)
        this->addZones(zoneSizes, zoneNames);
    }
    else
    {
        DynamicList<Face> dynFaces(faces.size());
        DynamicList<label> dynZones(std::move(zones));

        for (const face& f : faces)
        {
            dynFaces.append(Face(f));
        }

        // Count
        labelList zoneSizes(nZones, Zero);
        for (const label zonei : dynZones)
        {
            zoneSizes[zonei]++;
        }

        this->sortFacesAndStore(dynFaces, dynZones, dynElemId, sorted);

        // Add zones (retaining empty ones)
        this->addZones(zoneSizes, zoneNames);
    }
    this->addZonesToFaces(); // for labelledTri

    // transfer to normal lists
    this->storedPoints().transfer(reader.points());

    return true;
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    IOstreamOption,
    const dictionary& options
)
{
    const UList<point>& pointLst = surf.points();
    const UList<Face>&   faceLst = surf.surfFaces();
    const UList<label>&  faceMap = surf.faceMap();

    const surfZoneList zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    vtk::outputOptions opts = formatOptions(options);

    std::ofstream os(filename, std::ios::binary);

    autoPtr<vtk::formatter> format = opts.newFormatter(os);

    writeHeader(format(), pointLst);

    if (useFaceMap)
    {
        // connectivity count without additional storage (done internally)
        label nConnectivity = 0;
        for (const Face& f : faceLst)
        {
            nConnectivity += f.size();
        }

        vtk::legacy::beginPolys
        (
            format().os(),
            faceLst.size(),
            nConnectivity
        );

        label faceIndex = 0;
        for (const surfZone& zone : zones)
        {
            forAll(zone, i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                format().write(f.size());  // The size prefix
                vtk::writeList(format(), f);
            }
        }

        format().flush();
    }
    else
    {
        // Easy to write polys without a faceMap
        writePolys(format(), faceLst);
    }

    // Write regions (zones) as CellData
    if (zones.size() > 1)
    {
        writeCellData(format(), zones);
    }
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf,
    IOstreamOption,
    const dictionary& options
)
{
    vtk::outputOptions opts = formatOptions(options);

    std::ofstream os(filename, std::ios::binary);

    autoPtr<vtk::formatter> format = opts.newFormatter(os);

    writeHeader(format(), surf.points());

    // Easy to write polys without a faceMap
    writePolys(format(), surf.surfFaces());

    // Write regions (zones) as CellData
    writeCellData(format(), surf.zoneIds());
}


// ************************************************************************* //
