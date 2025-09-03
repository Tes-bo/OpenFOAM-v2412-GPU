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

// OpenFOAM includes
#include "functionObjectLine.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

// VTK Readers
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeName(functionObjectLine);
    addToRunTimeSelectionTable(pathline, functionObjectLine, dictionary);
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

    return dataset;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectLine::functionObjectLine
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    pathline(parent, dict, colours),
    functionObjectBase(parent, dict, colours),
    actor_()
{
    actor_ = vtkSmartPointer<vtkActor>::New();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectLine::~functionObjectLine()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::functionObjectLine::
addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    // Currently master-only
    if (!visible_ || !renderer || !Pstream::master())
    {
        return;
    }

    fileName fName = getFileName("file", fieldName_);
    if (fName.empty())
    {
        WarningInFunction
            << "Unable to read file name from function object "
            << functionObjectName_ << " for field " << fieldName_
            << ". Line will not be processed"
            << endl;
        return;
    }


    auto polyData = getPolyDataFile(fName);

    if (!polyData || polyData->GetNumberOfPoints() == 0)
    {
        WarningInFunction
            << "Could not read "<< fName << nl
            << "Only VTK (.vtp, .vtk) files are supported"
            << endl;
        return;
    }

    DebugInfo << "    Resolved lines " << fName << endl;


    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    setField
    (
        position,
        fieldName_,
        queryFieldAssociation(fieldName_, polyData),
        mapper,
        renderer
    );

    actor_->SetMapper(mapper);

    addLines(position, actor_, polyData);

    renderer->AddActor(actor_);
}


void Foam::functionObjects::runTimePostPro::functionObjectLine::updateActors
(
    const scalar position
)
{
    actor_->GetProperty()->SetLineWidth(2);
    actor_->GetProperty()->SetOpacity(opacity(position));
}


bool Foam::functionObjects::runTimePostPro::functionObjectLine::clear()
{
    if (functionObjectBase::clear())
    {
        return removeFile("file", fieldName_);
    }

    return false;
}


// ************************************************************************* //
