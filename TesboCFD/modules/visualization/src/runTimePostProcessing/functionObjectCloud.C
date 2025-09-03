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
#include "functionObjectCloud.H"
#include "fvMesh.H"
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
    defineTypeName(functionObjectCloud);
    addToRunTimeSelectionTable(pointData, functionObjectCloud, dictionary);
}
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

static vtkSmartPointer<vtkPolyData> getPolyDataFile(const Foam::fileName& fName)
{
    // Very simple - we only support vtp files, which are expected to have
    // the scaling and colouring fields.

    vtkSmartPointer<vtkPolyData> dataset;

    if (fName.ext() == "vtp")
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

Foam::functionObjects::runTimePostPro::functionObjectCloud::functionObjectCloud
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    pointData(parent, dict, colours),
    functionObjectBase(parent, dict, colours),
    cloudName_(dict.get<word>("cloud")),
    inputFileName_(),
    colourFieldName_(dict.get<word>("colourField")),
    actor_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectCloud::
~functionObjectCloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::functionObjectCloud::
addGeometryFromFile
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return false;
    }

    vtkSmartPointer<vtkPolyData> polyData;

    bool good = true;

    // The vtkCloud stores 'file' via the stateFunctionObject
    // (lookup by cloudName).
    // It only generates VTP format, which means there is a single file
    // containing all fields.

    if (Pstream::master())
    {
        inputFileName_ = getFileName("file", cloudName_);

        if (inputFileName_.size())
        {
            polyData = getPolyDataFile(inputFileName_);

            if (!polyData || polyData->GetNumberOfPoints() == 0)
            {
                good = false;

                WarningInFunction
                    << "Could not read "<< inputFileName_ << nl
                    << "Only VTK (.vtp) files are supported"
                    << endl;
            }
            else
            {
                DebugInfo
                    << "    Resolved cloud file " << inputFileName_ << endl;
            }
        }
        else
        {
            good = false;

            WarningInFunction
                << "Unable to find function object " << functionObjectName_
                << " output for field " << fieldName_
                << ". Cloud will not be processed"
                << endl;
        }
    }
    else
    {
        inputFileName_.clear();
    }

    reduce(good, andOp<bool>());

    if (!good)
    {
        return false;
    }

    // Only render on master
    if (!renderer || !Pstream::master())
    {
        return true;
    }


    // Rendering

    actor_ = vtkSmartPointer<vtkActor>::New();

    {
        fieldSummary scaleFieldInfo =
            queryFieldSummary(fieldName_, polyData);

        fieldSummary colourFieldInfo =
            queryFieldSummary(colourFieldName_, polyData);

        DebugInfo
            << "    Field " << fieldName_ << ' ' << scaleFieldInfo.info() << nl
            << "    Field " << colourFieldName_ << ' ' << colourFieldInfo.info()
            << endl;


        // No reduction

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        actor_->SetMapper(mapper);

        /// dataset->Print(std::cout);

        addGlyphs
        (
            position,
            fieldName_, scaleFieldInfo,         // scaling
            colourFieldName_, colourFieldInfo,  // colour
            maxGlyphLength_,
            polyData,
            actor_,
            renderer
        );

        renderer->AddActor(actor_);
    }

    return true;
}


void Foam::functionObjects::runTimePostPro::functionObjectCloud::
addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    // File source
    addGeometryFromFile(position, renderer);
}


void Foam::functionObjects::runTimePostPro::functionObjectCloud::updateActors
(
    const scalar position
)
{
    if (actor_)
    {
        const vector colour = pointColour_->value(position);

        vtkProperty* prop = actor_->GetProperty();

        prop->SetOpacity(opacity(position));

        prop->SetColor(colour[0], colour[1], colour[2]);
    }
}


bool Foam::functionObjects::runTimePostPro::functionObjectCloud::clear()
{
    if (functionObjectBase::clear())
    {
        if (inputFileName_.size() && Foam::rm(inputFileName_))
        {
            inputFileName_.clear();
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
