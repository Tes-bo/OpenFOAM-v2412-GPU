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

#include "triSurface.H"
#include "Fstream.H"
#include "Time.H"
#include "boundBox.H"
#include "bitSet.H"
#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "MeshedSurfaceProxy.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::wordHashSet Foam::triSurface::readTypes()
{
    wordHashSet known
    (
        UnsortedMeshedSurface<labelledTri>::readTypes()
      | MeshedSurface<labelledTri>::readTypes()
    );

    // Additional hard-coded entry points, but do not need
    // {"stl", "stlb"}
    // since they are shadowed by *MeshedSurface

    known.insert("ftr");

    return known;
}


Foam::wordHashSet Foam::triSurface::writeTypes()
{
    wordHashSet known
    (
        MeshedSurfaceProxy<labelledTri>::writeTypes()
    );

    // Additional hard-coded entry points, but do not need
    // {"gts", "stl", "stlb"}
    // since they are shadowed by MeshedSurfaceProxy

    known.insert("ftr");

    return known;
}


bool Foam::triSurface::canReadType(const word& fileType, bool verbose)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        readTypes(),
        fileType,
        verbose,
        "reading"
    );
}


bool Foam::triSurface::canWriteType(const word& fileType, bool verbose)
{
    return fileFormats::surfaceFormatsCore::checkSupport
    (
        writeTypes(),
        fileType,
        verbose,
        "writing"
    );
}


bool Foam::triSurface::canRead(const fileName& name, bool verbose)
{
    const word ext =
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );

    return canReadType(ext, verbose);
}


Foam::fileName Foam::triSurface::relativeFilePath
(
    const IOobject& io,
    const fileName& f,
    const bool isGlobal
)
{
    return fileFormats::surfaceFormatsCore::relativeFilePath(io, f, isGlobal);
}


Foam::fileName Foam::triSurface::checkFile
(
    const IOobject& io,
    const bool isGlobal
)
{
    return fileFormats::surfaceFormatsCore::checkFile(io, isGlobal);
}


Foam::fileName Foam::triSurface::checkFile
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
{
    return fileFormats::surfaceFormatsCore::checkFile(io, dict, isGlobal);
}


Foam::fileName Foam::triSurface::findFile
(
    const IOobject& io,
    const bool isGlobal
)
{
    return fileFormats::surfaceFormatsCore::findFile(io, isGlobal);
}


Foam::fileName Foam::triSurface::findFile
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
{
    return fileFormats::surfaceFormatsCore::findFile(io, dict, isGlobal);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triSurface::readNative(Istream& is)
{
    // Read triangles, points from Istream
    is  >> patches_ >> storedPoints() >> storedFaces();

    return true;
}


void Foam::triSurface::writeNative(Ostream& os) const
{
    os  << patches() << nl;

    //Note: Write with global point numbering
    os  << points() << nl
        << static_cast<const List<labelledTri>&>(*this) << nl;

    os.check(FUNCTION_NAME);
}


bool Foam::triSurface::read
(
    const fileName& name,
    const word& fileType,
    const bool check
)
{
    if (check && !exists(name))
    {
        FatalErrorInFunction
            << "No such file " << name << nl
            << exit(FatalError);
    }

    this->clear();
    transfer(*New(name, fileType));
    return true;
}


void Foam::triSurface::write
(
    const fileName& name,
    const word& fileType,
    const bool sortByRegion
) const
{
    if (fileType.empty())
    {
        // Handle empty/missing type

        const word ext(name.ext());

        if (ext.empty())
        {
            FatalErrorInFunction
                << "Cannot determine format from filename" << nl
                << "    " << name << nl
                << exit(FatalError);
        }

        write(name, ext, sortByRegion);
        return;
    }


    // Hard-coded writers

    if (fileType == "ftr")
    {
        OFstream os(name);
        writeNative(os);
    }
    else if (fileType == "stl")
    {
        writeSTLASCII(name, sortByRegion);
    }
    else if (fileType == "stlb")
    {
        writeSTLBINARY(name);
    }
    else if (fileType == "gts")
    {
        writeGTS(name, sortByRegion);
    }
    else if (MeshedSurfaceProxy<labelledTri>::canWriteType(fileType))
    {
        labelList faceMap;
        List<surfZone> zoneLst = this->sortedZones(faceMap);

        MeshedSurfaceProxy<labelledTri> proxy
        (
            this->points(),
            this->surfFaces(),
            zoneLst,
            faceMap
        );

        proxy.write(name, fileType);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown surface format " << fileType
            << " for writing file " << name << nl
            << "Valid types:" << nl
            << "    " << flatOutput(writeTypes().sortedToc()) << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurface::triSurface(Istream& is)
:
    triSurface()
{
    readNative(is);

    setDefaultPatches();
}


Foam::triSurface::triSurface(const Time& d)
:
    triSurface()
{
    IFstream is
    (
        d.path()/triSurfInstance(d)/typeName/(d.caseName() + ".ftr")
    );

    readNative(is);

    setDefaultPatches();
}


Foam::triSurface::triSurface
(
    const IOobject& io,
    const dictionary& dict,
    const bool isGlobal
)
:
    triSurface()
{
    fileName fName(checkFile(io, dict, isGlobal));

    read(fName, dict.getOrDefault<word>("fileType", word::null));

    scalePoints(dict.getOrDefault<scalar>("scale", 0));

    setDefaultPatches();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::write
(
    const fileName& name,
    const bool sortByRegion
) const
{
    write(name, name.ext(), sortByRegion);
}


void Foam::triSurface::write(Ostream& os) const
{
    writeNative(os);
}


void Foam::triSurface::write(const Time& d) const
{
    OFstream os
    (
        d.path()/triSurfInstance(d)/typeName/(d.caseName() + ".ftr")
    );

    writeNative(os);
}


void Foam::triSurface::writeStats(Ostream& os) const
{
    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves.

    bitSet pointIsUsed(points().size());

    boundBox bb;
    labelHashSet regionsUsed;

    for (const auto& f : *this)
    {
        regionsUsed.insert(f.region());

        for (const label pointi : f)
        {
            if (pointIsUsed.set(pointi))
            {
                bb.add(points()[pointi]);
            }
        }
    }

    os  << "Triangles    : " << size()
        << " in " << regionsUsed.size() <<  " region(s)" << nl
        << "Vertices     : " << pointIsUsed.count() << nl
        << "Bounding Box : " << bb << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, triSurface& s)
{
    s.clearOut();
    s.readNative(is);
    s.setDefaultPatches();
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const triSurface& s)
{
    s.writeNative(os);
    return os;
}


// ************************************************************************* //
