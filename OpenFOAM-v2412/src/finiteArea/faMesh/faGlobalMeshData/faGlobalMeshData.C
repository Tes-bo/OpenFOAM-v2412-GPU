/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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

Author
    Hrvoje Jasak

\*----------------------------------------------------------------------------*/

#include "faGlobalMeshData.H"
#include "faMesh.H"
#include "globalMeshData.H"
#include "processorFaPatch.H"
#include "processorTopologyNew.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faGlobalMeshData::faGlobalMeshData(const faMesh& mesh)
:
    mesh_(mesh),
    processorTopology_
    (
        processorTopology::New<processorFaPatch>
        (
            mesh.boundary(),
            UPstream::worldComm
        )
    ),
    nGlobalPoints_(-1),
    sharedPointLabels_(),
    sharedPointAddr_()
{
    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// A non-default destructor since we had incomplete types in the header
Foam::faGlobalMeshData::~faGlobalMeshData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faMesh& Foam::faGlobalMeshData::mesh() const noexcept
{
    return mesh_;
}


void Foam::faGlobalMeshData::updateMesh()
{
    // The polyMesh globalMeshData
    const globalMeshData& pData = mesh_.mesh().globalData();

    const labelList& polyMeshSharedPointLabels =
        pData.sharedPointLabels();

    const labelList& polyMeshSharedPointAddr =
        pData.sharedPointAddr();

    labelHashSet sharedPointLabels;

    labelField globalList(pData.nGlobalPoints(), Zero);

    forAll(mesh_.boundary(), patchI)
    {
        const faPatch& fap = mesh_.boundary()[patchI];

        if (isA<processorFaPatch>(fap))
        {
            const labelList& localPointLabels = fap.pointLabels();

            forAll(localPointLabels, pointI)
            {
                label polyMeshPoint =
                    mesh_.patch().meshPoints()[localPointLabels[pointI]];

                const label sharedPolyMeshPoint =
                    polyMeshSharedPointLabels.find(polyMeshPoint);

                if
                (
                    sharedPolyMeshPoint != -1
                 && !sharedPointLabels.found(localPointLabels[pointI])
                )
                {
                    globalList[polyMeshSharedPointAddr[sharedPolyMeshPoint]]
                        += 1;

                    sharedPointLabels.insert(localPointLabels[pointI]);
                }
            }
        }
    }

    sharedPointLabels_ = sharedPointLabels.toc();

    Pstream::combineReduce(globalList, plusEqOp<labelField>());

    nGlobalPoints_ = 0;
    for (label i=0; i<globalList.size(); ++i)
    {
        if (globalList[i] > 0)
        {
            globalList[i] = ++nGlobalPoints_;
        }
    }

    sharedPointAddr_.setSize(sharedPointLabels_.size());
    forAll(sharedPointAddr_, pointI)
    {
        const label polyMeshSharedPointIndex =
            polyMeshSharedPointLabels.find
            (
                mesh_.patch().meshPoints()[sharedPointLabels_[pointI]]
            );

        sharedPointAddr_[pointI] =
            globalList[polyMeshSharedPointAddr[polyMeshSharedPointIndex]]
          - 1;
    }
}


// ************************************************************************* //
