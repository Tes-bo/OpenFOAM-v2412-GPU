/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "columnAverage.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "meshStructure.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(columnAverage, 0);
    addToRunTimeSelectionTable(functionObject, columnAverage, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::meshStructure&
Foam::functionObjects::columnAverage::meshAddressing(const polyMesh& mesh) const
{
    if (!meshStructurePtr_)
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        // Count
        label sz = 0;
        for (const label patchi : patchIDs_)
        {
            sz += pbm[patchi].size();
        }

        // Fill
        labelList meshFaces(sz);
        sz = 0;
        for (const label patchi : patchIDs_)
        {
            label start = pbm[patchi].start();
            label size = pbm[patchi].size();
            for (label i = 0; i < size; ++i)
            {
                meshFaces[sz++] = start+i;
            }
        }

        if (sz == 0)
        {
            // TODO: If source patch is a cyclic it may have have been
            // converted to a processorCyclic for parallel runs

            WarningInFunction
                << "Requested patches have zero faces"
                << endl;
        }

        uindirectPrimitivePatch uip
        (
            UIndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        );

        globalFaces_.reset(new globalIndex(uip.size()));
        globalEdges_.reset(new globalIndex(uip.nEdges()));
        globalPoints_.reset(new globalIndex(uip.nPoints()));
        meshStructurePtr_.reset
        (
            new meshStructure
            (
                mesh,
                uip,
                globalFaces_(),
                globalEdges_(),
                globalPoints_()
            )
        );
    }

    return *meshStructurePtr_;
}


const Foam::word Foam::functionObjects::columnAverage::averageName
(
    const word& fieldName
) const
{
    return scopedName("columnAverage(" + fieldName + ")");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::columnAverage::columnAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchIDs_(),
    fieldSet_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::columnAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    patchIDs_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc();

    fieldSet_.read(dict);

    return true;
}


bool Foam::functionObjects::columnAverage::execute()
{
    // Make fields up to date with current selection
    fieldSet_.updateSelection();

    for (const word& fieldName : fieldSet_.selectionNames())
    {
        columnAverageField<scalar>(fieldName);
        columnAverageField<vector>(fieldName);
        columnAverageField<sphericalTensor>(fieldName);
        columnAverageField<symmTensor>(fieldName);
        columnAverageField<tensor>(fieldName);
    }

    return true;
}


bool Foam::functionObjects::columnAverage::write()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const regIOobject* obj =
            obr_.cfindObject<regIOobject>(averageName(fieldName));

        if (obj)
        {
            obj->write();
        }
    }

    return true;
}


// ************************************************************************* //
