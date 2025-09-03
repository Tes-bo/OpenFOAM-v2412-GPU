/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Face>
Foam::autoPtr<Foam::MeshedSurface<Face>>
Foam::MeshedSurface<Face>::New
(
    const fileName& name,
    const word& fileType,
    bool mandatory
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

        return New(name, ext, mandatory);
    }
    else if (fileType == "gz")
    {
        // Degenerate call
        return New(name.lessExt(), name.stem().ext(), mandatory);
    }
    else if (ext == "gz")
    {
        // Handle trailing "gz" on file name
        return New(name.lessExt(), fileType, mandatory);
    }

    // if (check && !exists(name))
    // {
    //     FatalErrorInFunction
    //         << "No such file " << name << nl
    //         << exit(FatalError);
    // }


    DebugInFunction
        << "Construct MeshedSurface (" << fileType << ")\n";

    auto* ctorPtr = fileExtensionConstructorTable(fileType);

    if (ctorPtr)
    {
        return autoPtr<MeshedSurface<Face>>(ctorPtr(name));
    }


    // Delegate to friend if possible
    const wordHashSet delegate(FriendType::readTypes());

    if (delegate.found(fileType))
    {
        // OK, can create indirectly
        auto surf = autoPtr<MeshedSurface<Face>>::New();
        surf->transfer(*FriendType::New(name, fileType));

        return surf;
    }
    else if (mandatory)
    {
        FatalErrorInFunction
            << "Unknown surface format " << fileType << nl << nl
            << "Valid types:" << nl
            << flatOutput((delegate | readTypes()).sortedToc()) << nl
            << exit(FatalError);
    }

    // Failed, but was optional
    return nullptr;
}


template<class Face>
Foam::autoPtr<Foam::MeshedSurface<Face>>
Foam::MeshedSurface<Face>::New(const fileName& name)
{
    if (name.has_ext("gz"))
    {
        // Handle trailing "gz" on file name
        return New(name.lessExt(), name.stem().ext());
    }

    return New(name, name.ext());
}


// ************************************************************************* //
