/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

// OpenFOAM includes
#include "geometrySurface.H"
#include "stringOps.H"
#include "foamVtkTools.H"
#include "MeshedSurfaces.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

// VTK Readers
#include "vtkOBJReader.h"
#include "vtkSTLReader.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeNameAndDebug(geometrySurface, 0);
    addToRunTimeSelectionTable(surface, geometrySurface, dictionary);
}
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

static vtkSmartPointer<vtkPolyData> getPolyDataFile(const Foam::fileName& fName)
{
    // Not extremely elegant...
    vtkSmartPointer<vtkPolyData> dataset;

    if ("vtk" == fName.ext())
    {
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    if ("vtp" == fName.ext())
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    if ("obj" == fName.ext())
    {
        auto reader = vtkSmartPointer<vtkOBJReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    if ("stl" == fName.ext() || "stlb" == fName.ext())
    {
        auto reader = vtkSmartPointer<vtkSTLReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }


    // Fallback to using OpenFOAM to read the surface and convert afterwards
    Foam::meshedSurface surf(fName);

    dataset = Foam::vtk::Tools::Patch::mesh(surf);

    dataset->GetCellData()->SetNormals
    (
        Foam::vtk::Tools::Patch::faceNormals(surf)
    );

    return dataset;
}

} // End anonymous namespace


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::geometrySurface::addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer,
    const fileName& fName
) const
{
    // Master-only, since that is where the files are.
    if (!visible_ || !renderer || !Pstream::master())
    {
        return;
    }

    if (representation_ == rtGlyph)
    {
        FatalErrorInFunction
            << "Glyph representation not available for " << typeName
            << " object" << exit(FatalError);
    }

    DebugInfo << "    Add geometry surface: " << fName << nl;

    auto surf = getPolyDataFile(fName);

    if (!surf || !surf->GetNumberOfPoints())
    {
        FatalErrorInFunction
            << "Could not read "<< fName << nl
            << exit(FatalError);
    }

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    mapper->ScalarVisibilityOff();

    mapper->SetInputData(surf);

    addFeatureEdges(renderer, surf);

    surfaceActor_->SetMapper(mapper);

    setRepresentation(surfaceActor_);

    renderer->AddActor(surfaceActor_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometrySurface::geometrySurface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    surface(parent, dict, colours),
    fileNames_()
{
    dict.readEntry("files", fileNames_);
}


Foam::functionObjects::runTimePostPro::geometrySurface::geometrySurface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours,
    const List<fileName>& fileNames
)
:
    surface(parent, dict, colours),
    fileNames_(fileNames)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::geometrySurface::addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return;
    }

    for (fileName f : fileNames_)  // Use a copy
    {
        f.expand();
        addGeometryToScene(position, renderer, f);
    }
}


void Foam::functionObjects::runTimePostPro::geometrySurface::updateActors
(
    const scalar position
)
{
    if (!visible_)
    {
        return;
    }

    surface::updateActors(position);
}


bool Foam::functionObjects::runTimePostPro::geometrySurface::clear()
{
    // Note: do not remove geometry files
    // - often static files used for other purposes as well (eg meshing)

    return true;
}


// ************************************************************************* //
