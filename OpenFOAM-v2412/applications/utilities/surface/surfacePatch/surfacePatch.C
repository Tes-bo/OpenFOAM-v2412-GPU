/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    surfacePatch

Group
    grpSurfaceUtilities

Description
    Patches (regionises) a surface using a user-selectable method.

    The current methods are either based on a geometric feature angle
    (a replacement for the surfaceAutoPatch functionality) or on intersecting
    a set of searchableSurfaces - this will re-triangulate the intersecting
    triangles and regonise the resulting triangles according to being
    inside/outside the surface.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "searchableSurfaces.H"
#include "searchableSurfaceModifier.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Add patches (regions) to a surface with a user-selectable method"
    );
    argList::noParallel();
    argList::addOption("dict", "file", "Alternative surfacePatchDict");

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("surfacePatchDict");
    #include "setSystemRunTimeDictionaryIO.H"
    const IOdictionary meshDict(dictIO);


    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            runTime.constant(),         // instance
            triSurfaceMesh::meshSubDir, // local
            runTime,                    // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshDict.subDict("geometry"),
        meshDict.getOrDefault("singleRegionName", true)
    );


    const dictionary& surfacesDict = meshDict.subDict("surfaces");

    for (const entry& dEntry : surfacesDict)
    {
        if (dEntry.isDict())
        {
            const word& surfName = dEntry.keyword();
            const dictionary& surfDict = dEntry.dict();

            // Look up surface
            searchableSurface& surf = allGeometry[surfName];


            bool changed = false;

            // Check for modifier on global level
            if (surfDict.found("type"))
            {
                autoPtr<searchableSurfaceModifier> modifier
                (
                    searchableSurfaceModifier::New
                    (
                        surfDict.get<word>("type"),
                        allGeometry,
                        surfDict
                    )
                );

                if (modifier().modify(identity(surf.regions().size()), surf))
                {
                    changed = true;
                }
            }

            // Check for modifier on a per-region level
            if (surfDict.found("regions"))
            {
                const dictionary& regionsDict = surfDict.subDict("regions");

                for (const entry& e : regionsDict)
                {
                    const wordRe regionName(e.keyword());
                    const dictionary& regionDict = e.dict();

                    labelList regionIDs
                    (
                        wordRes::matching(regionName, surf.regions())
                    );

                    autoPtr<searchableSurfaceModifier> modifier
                    (
                        searchableSurfaceModifier::New
                        (
                            regionDict.get<word>("type"),
                            allGeometry,
                            regionDict
                        )
                    );

                    if (modifier().modify(regionIDs, surf))
                    {
                        changed = true;
                    }
                }
            }


            if (changed)
            {
                const word oldName(surf.name());
                surf.rename(oldName.lessExt() + "_patched." + oldName.ext());

                Info<< "Writing repatched surface " << oldName << " to "
                    << surf.name() << nl << endl;

                surf.write();
            }
            else
            {
                Info<< "Surface " << surf.name() << " unchanged" << nl << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
