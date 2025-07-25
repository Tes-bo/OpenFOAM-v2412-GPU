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

#include "triSurface.H"
#include "Fstream.H"
#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::triSurface>
Foam::triSurface::New
(
    const fileName& name,
    const word& fileType
)
{
    const word ext(name.ext());

    if (fileType.empty())
    {
        // Handle empty/missing type

        if (ext.empty())
        {
            FatalErrorInFunction
                << "Cannot determine format from filename" << nl
                << "    " << name << nl
                << exit(FatalError);
        }

        return New(name, ext);
    }
    else if (fileType == "gz")
    {
        // Degenerate call
        return New(name.lessExt(), name.stem().ext());
    }
    else if (ext == "gz")
    {
        // Handle trailing "gz" on file name
        return New(name.lessExt(), fileType);
    }

    // if (check && !exists(name))
    // {
    //     FatalErrorInFunction
    //         << "No such file " << name << nl
    //         << exit(FatalError);
    // }


    // Hard-coded readers
    if (fileType == "ftr")
    {
        auto surf = autoPtr<triSurface>::New();

        IFstream is(name);
        surf->readNative(is);
        return surf;
    }
    else if (fileType == "stl")
    {
        auto surf = autoPtr<triSurface>::New();

        surf->readSTL(name);  // ASCII
        return surf;
    }
    else if (fileType == "stlb")
    {
        auto surf = autoPtr<triSurface>::New();

        surf->readSTL(name, true); // Force BINARY
        return surf;
    }

    {
        // UnsortedMeshedSurface
        using proxyType = UnsortedMeshedSurface<labelledTri>;
        if (proxyType::readTypes().found(fileType))
        {
            auto surf = autoPtr<triSurface>::New();

            surf->transfer(*proxyType::New(name, fileType));
            return surf;
        }
    }

    // MeshedSurface
    {
        using proxyType = MeshedSurface<labelledTri>;
        if (proxyType::readTypes().found(fileType))
        {
            auto surf = autoPtr<triSurface>::New();

            surf->transfer(*proxyType::New(name, fileType));
            return surf;
        }
    }

    {
        FatalErrorInFunction
            << "Unknown surface format " << fileType
            << " for reading file " << name << nl
            << "Valid types:" << nl
            << "    " << flatOutput(readTypes().sortedToc()) << nl
            << exit(FatalError);
    }

    // Failed
    return nullptr;
}


Foam::autoPtr<Foam::triSurface>
Foam::triSurface::New(const fileName& name)
{
    if (name.has_ext("gz"))
    {
        // Handle trailing "gz" on file name
        return New(name.lessExt(), name.stem().ext());
    }

    return New(name, name.ext());
}


// ************************************************************************* //
