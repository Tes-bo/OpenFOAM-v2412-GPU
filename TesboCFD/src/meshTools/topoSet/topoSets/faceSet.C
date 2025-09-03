/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "faceSet.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "mapDistributePolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(faceSet);
    addToRunTimeSelectionTable(topoSet, faceSet, word);
    addToRunTimeSelectionTable(topoSet, faceSet, size);
    addToRunTimeSelectionTable(topoSet, faceSet, set);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceSet::faceSet(const IOobject& io)
:
    topoSet(io, typeName)
{}


Foam::faceSet::faceSet(const IOobject& io, const Foam::zero)
:
    topoSet(io, Foam::zero{})
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt,
    IOobjectOption::registerOption reg
)
:
    topoSet(mesh, typeName, name, rOpt, wOpt, reg)
{
    check(mesh.nFaces());  // Valid range?
}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const label initialCapacity,
    IOobjectOption::writeOption wOpt
)
:
    topoSet(mesh, name, initialCapacity, wOpt)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    IOobjectOption::writeOption wOpt
)
:
    topoSet(mesh, name, set, wOpt)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& labels,
    IOobjectOption::writeOption wOpt
)
:
    topoSet(mesh, name, labels, wOpt)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    labelHashSet&& labels,
    IOobjectOption::writeOption wOpt
)
:
    topoSet(mesh, name, std::move(labels), wOpt)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const labelUList& labels,
    IOobjectOption::writeOption wOpt
)
:
    topoSet(mesh, name, labels, wOpt)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::labelHashSet Foam::faceSet::readContents
(
    const polyMesh& mesh,
    const word& name
)
{
    faceSet reader
    (
        topoSet::findIOobject(mesh, name, IOobjectOption::NO_REGISTER),
        Foam::zero{}
    );

    labelHashSet labels;
    reader.readIOcontents(typeName, labels);
    reader.checkLabels(labels, mesh.nFaces());  // Valid range?

    return labels;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceSet::sync(const polyMesh& mesh)
{
    labelHashSet& labels = *this;

    // Convert to boolList
    // TBD: could change to using bitSet for the synchronization

    const label len = mesh.nFaces();

    boolList contents(len, false);

    for (const label facei : labels)
    {
        contents.set(facei);
    }

    syncTools::syncFaceList(mesh, contents, orEqOp<bool>());


    // Update labelHashSet

    for (label i=0; i < len; ++i)
    {
        if (contents.test(i))
        {
            labels.set(i);
        }
    }
}


Foam::label Foam::faceSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


void Foam::faceSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseFaceMap());
}


void Foam::faceSet::distribute(const mapDistributePolyMesh& map)
{
    labelHashSet& labels = *this;

    boolList contents(map.nOldFaces(), false);

    for (const label facei : labels)
    {
        contents.set(facei);
    }

    map.distributeFaceData(contents);

    // The new length
    const label len = contents.size();

    // Count - as per BitOps::count(contents)
    label n = 0;
    for (label i=0; i < len; ++i)
    {
        if (contents.test(i))
        {
            ++n;
        }
    }

    // Update labelHashSet

    labels.clear();
    labels.reserve(n);

    for (label i=0; i < len; ++i)
    {
        if (contents.test(i))
        {
            labels.set(i);
        }
    }
}


void Foam::faceSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.faceCentres(), maxLen);
}


// ************************************************************************* //
