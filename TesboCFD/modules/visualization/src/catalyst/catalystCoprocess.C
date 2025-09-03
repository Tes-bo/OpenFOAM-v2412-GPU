/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
    Copyright (C) 2018 CINECA
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

#include "catalystCoprocess.H"
#include "Time.H"

// VTK includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>

// Catalyst includes
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{
    defineTypeNameAndDebug(coprocess, 0);
}
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class DataType>
bool Foam::catalyst::coprocess::processImpl
(
    const dataQuery& dataq,
    HashTable<vtkSmartPointer<DataType>>& outputs
)
{
    auto* descrip = dataq.description();

    if (!coproc_->RequestDataDescription(descrip))
    {
        return false;
    }

    for (const word& channel : dataq.names())
    {
        if (outputs.contains(channel))
        {
            auto* input = descrip->GetInputDescriptionByName(channel.c_str());

            if (input && input->GetIfGridIsNecessary())
            {
                input->SetGrid(outputs[channel]);
            }
        }
    }

    coproc_->CoProcess(descrip);
    return true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::catalyst::coprocess::~coprocess()
{
    // As per stop(), but without output
    if (coproc_)
    {
        coproc_->Delete();
        coproc_ = nullptr;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::catalyst::coprocess::good() const
{
    return coproc_;
}


void Foam::catalyst::coprocess::stop()
{
    if (coproc_)
    {
        coproc_->Delete();
        coproc_ = nullptr;
        Info<< "Stop catalyst" << endl;
    }
}


void Foam::catalyst::coprocess::reset(const fileName& outputDir)
{
    if (coproc_ == nullptr)
    {
        coproc_ = vtkCPProcessor::New();
        coproc_->Initialize(outputDir.c_str());
        Info<< "Connecting ParaView Catalyst..." << endl;
    }
    else
    {
        coproc_->RemoveAllPipelines();

        if (outputDir == coproc_->GetWorkingDirectory())
        {
            Info<< "Rebinding ParaView Catalyst..." << endl;
        }
        else
        {
            // Changed working directory ... redo everything.
            coproc_->Delete();
            coproc_ = nullptr;

            reset(outputDir);
        }
    }
}


void Foam::catalyst::coprocess::reset
(
    const fileName& outputDir,
    const UList<string>& scripts
)
{
    reset(outputDir);

    int nscript = 0;
    for (const auto& script : scripts)
    {
        Info<< "Adding pipeline[" << nscript << "] " << script << endl;
        ++nscript;

        vtkNew<vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(script.c_str());
        coproc_->AddPipeline(pipeline.GetPointer());
    }

    // Do something different with (nscript == 0) ?
}


Foam::label Foam::catalyst::coprocess::query(dataQuery& dataq)
{
    return dataq.query(coproc_);
}


bool Foam::catalyst::coprocess::process
(
    const dataQuery& dataq,
    outputChannels& outputs
)
{
    return processImpl(dataq, outputs);
}


// ************************************************************************* //
