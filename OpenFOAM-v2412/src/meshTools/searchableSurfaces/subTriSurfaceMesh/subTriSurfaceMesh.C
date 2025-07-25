/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "subTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(subTriSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        subTriSurfaceMesh,
        dict
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::subTriSurfaceMesh::patchNames(const triSurface& s)
{
    const auto& patches = s.patches();

    wordList names(patches.size());
    forAll(patches, patchi)
    {
        names[patchi] = patches[patchi].name();
    }
    return names;
}


Foam::labelList Foam::subTriSurfaceMesh::selectedRegions
(
    const triSurface& s,
    const wordRes& regionNameMatcher
)
{
    const wordList names(patchNames(s));

    labelList regionIds(names.size());

    label count = 0;

    forAll(names, regioni)
    {
        if (regionNameMatcher.match(names[regioni]))
        {
            regionIds[count++] = regioni;
        }
    }
    regionIds.setSize(count);

    return regionIds;
}


Foam::triSurface Foam::subTriSurfaceMesh::subset
(
    const IOobject& io,
    const dictionary& dict
)
{
    const word subGeomName(dict.get<word>("surface"));

    const triSurfaceMesh& s =
        io.db().lookupObject<triSurfaceMesh>(subGeomName);

    const wordRes regionNames(dict.get<wordRes>("patches"));

    labelList regionMap(selectedRegions(s, regionNames));

    if (regionMap.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Found no regions in triSurface matching " << regionNames
            << ". Valid regions are " << patchNames(s)
            << exit(FatalIOError);
    }

    labelList reverseRegionMap(s.patches().size(), -1);
    forAll(regionMap, i)
    {
        reverseRegionMap[regionMap[i]] = i;
    }

    boolList isSelected(s.size(), false);
    forAll(s, triI)
    {
        if (reverseRegionMap[s.triSurface::operator[](triI).region()] != -1)
        {
            isSelected[triI] = true;
        }
    }

    return s.subsetMesh(isSelected);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subTriSurfaceMesh::subTriSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurfaceMesh(io, subset(io, dict))
{}


// ************************************************************************* //
