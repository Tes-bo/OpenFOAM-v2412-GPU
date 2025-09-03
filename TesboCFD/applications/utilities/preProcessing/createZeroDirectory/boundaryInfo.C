/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "boundaryInfo.H"
#include "Time.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOPtrList<Foam::entry> Foam::boundaryInfo::readBoundaryDict
(
    const Time& runTime,
    const word& regionName
) const
{
    Info<< "    Reading mesh boundaries" << endl;

    const_cast<word&>(IOPtrList<entry>::typeName) = polyBoundaryMesh::typeName;
    IOPtrList<entry> boundaryPatchList
    (
        IOobject
        (
            "boundary",
            runTime.findInstance(polyMesh::meshDir(regionName), "boundary"),
            polyMesh::meshDir(regionName),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    // remove zero-sized patches
    PtrList<entry> boundaryPatchListNew;
    forAll(boundaryPatchList, patchI)
    {
        const dictionary& dict = boundaryPatchList[patchI].dict();
        const word pType = dict.get<word>("type");
        bool procPatch = pType == processorPolyPatch::typeName;

        bool addPatch = true;
        if (!procPatch)
        {
            label nFaces = dict.get<label>("nFaces");
            if (returnReduceAnd(nFaces == 0))
            {
                addPatch = false;
            }
        }

        if (addPatch)
        {
            boundaryPatchListNew.append(boundaryPatchList[patchI].clone());
        }
    }

    boundaryPatchList.transfer(boundaryPatchListNew);

    return boundaryPatchList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryInfo::boundaryInfo(const Time& runTime, const word& regionName)
:
    boundaryDict_(readBoundaryDict(runTime, regionName)),
    names_(),
    types_(),
    constraint_(),
    groups_(),
    allGroupNames_()
{
    names_.setSize(boundaryDict_.size());
    types_.setSize(boundaryDict_.size());
    constraint_.setSize(boundaryDict_.size(), false);
    groups_.setSize(boundaryDict_.size());

    forAll(boundaryDict_, patchI)
    {
        const dictionary& dict = boundaryDict_[patchI].dict();

        names_[patchI] = dict.dictName();
        dict.readEntry("type", types_[patchI]);
        if (polyPatch::constraintType(types_[patchI]))
        {
            constraint_[patchI] = true;
        }

        if (dict.readIfPresent("inGroups", groups_[patchI]))
        {
            allGroupNames_.insert(groups_[patchI]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::boundaryInfo::names() const
{
    return names_;
}


const Foam::wordList& Foam::boundaryInfo::types() const
{
    return types_;
}


const Foam::boolList& Foam::boundaryInfo::constraint() const
{
    return constraint_;
}


const Foam::List<Foam::wordList>& Foam::boundaryInfo::groups() const
{
    return groups_;
}


const Foam::wordHashSet& Foam::boundaryInfo::allGroupNames() const
{
    return allGroupNames_;
}


void Foam::boundaryInfo::setType(const label patchI, const word& condition)
{
    if (constraint_[patchI])
    {
        // not overriding constraint types
        return;
    }

    if (regExp(".*[Mm]apped.*").match(types_[patchI]))
    {
        // ugly hack to avoid overriding mapped types
        return;
    }

    if (condition == "wall")
    {
        types_[patchI] = condition;
    }
    else
    {
        types_[patchI] = "patch";
    }

    dictionary& patchDict = boundaryDict_[patchI].dict();
    patchDict.add("type", types_[patchI], true);
}


void Foam::boundaryInfo::write() const
{
    boundaryDict_.write();
}


// ************************************************************************* //
