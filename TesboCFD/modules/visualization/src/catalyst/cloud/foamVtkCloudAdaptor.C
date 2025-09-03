/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "foamVtkCloudAdaptor.H"
#include "Cloud.H"
#include "IOField.H"
#include "predicates.H"

// Only used locally
#include "foamVtkCloudAdaptorTemplates.C"

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkMultiPieceDataSet.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{
    defineTypeNameAndDebug(cloudAdaptor, 0);
}
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class UnaryMatchPredicate>
vtkSmartPointer<vtkMultiPieceDataSet>
Foam::vtk::cloudAdaptor::getCloudImpl
(
    const objectRegistry& mesh,
    const word& cloudName,
    const UnaryMatchPredicate& matcher
)
{
    // All individual datasets are vtkMultiPieceDataSet for improved
    // handling downstream.

    label rank = 0;
    label nproc = 1;

    if (UPstream::parRun())
    {
        rank  = UPstream::myProcNo();
        nproc = UPstream::nProcs();
    }

    auto output = vtkSmartPointer<vtkMultiPieceDataSet>::New();

    const auto* objPtr = mesh.cfindObject<cloud>(cloudName);
    if (!objPtr)
    {
        return output;
    }

    objectRegistry obrTmp
    (
        IOobject
        (
            "vtk::catalyst::" + cloudName,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    objPtr->writeObjects(obrTmp);

    const auto* pointsPtr = obrTmp.cfindObject<vectorField>("position");

    vtkSmartPointer<vtkPolyData> vtkmesh;
    if (pointsPtr)
    {
        vtkmesh = vtk::Tools::Vertices(*pointsPtr);

        // Prevent any possible conversion of positions as a field
        obrTmp.filterKeys
        (
            [](const word& k)
            {
                return
                (
                    k.starts_with("position")
                 || k.starts_with("coordinate")
                );
            },
            true  // prune
        );

        // Restrict to specified fields
        obrTmp.filterKeys(matcher);

        convertLagrangianFields<label>(vtkmesh, obrTmp);
        convertLagrangianFields<scalar>(vtkmesh, obrTmp);
        convertLagrangianFields<vector>(vtkmesh, obrTmp);
    }
    else
    {
        // This should never occur
        vtkmesh = vtkSmartPointer<vtkPolyData>::New();
    }

    output->SetNumberOfPieces(nproc);
    output->SetPiece(rank, vtkmesh);

    return output;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::cloudAdaptor::cloudAdaptor(const fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkSmartPointer<vtkMultiPieceDataSet>
Foam::vtk::cloudAdaptor::getCloud
(
    const word& cloudName
) const
{
    return getCloudImpl(mesh_, cloudName, predicates::always());
}


vtkSmartPointer<vtkMultiPieceDataSet>
Foam::vtk::cloudAdaptor::getCloud
(
    const word& cloudName,
    const wordRes& selectFields
) const
{
    if (selectFields.empty())
    {
        return getCloud(cloudName);
    }

    return getCloudImpl(mesh_, cloudName, selectFields);
}


// ************************************************************************* //
