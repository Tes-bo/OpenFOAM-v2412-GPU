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
#include "runTimePostProcessing.H"
#include "geometrySurface.H"
#include "dictionary.H"
#include "pointData.H"
#include "pathline.H"
#include "surface.H"
#include "text.H"
#include "Time.H"
#include "sigFpe.H"
#include "polySurfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"
#include "vtkLight.h"

#include "vtkDummyController.h"

#ifdef FOAM_USING_VTK_MPI
# include "vtkMPICommunicator.h"
# include "vtkMPIController.h"
# include "vtkCompositedSynchronizedRenderers.h"
# include "vtkSynchronizedRenderWindows.h"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(runTimePostProcessing, 0);

    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        runTimePostProcessing,
        dictionary
    );

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static void addGeometryToScene
(
    PtrList<Type>& objects,
    const scalar position,
    vtkRenderer* renderer
)
{
    for (Type& obj : objects)
    {
        if (Pstream::master() || obj.parallel())
        {
            obj.addGeometryToScene(position, renderer);
        }
        else
        {
            obj.addGeometryToScene(position, nullptr);
        }
    }
}


template<class Type>
static void updateActors(PtrList<Type>& objects, const scalar position)
{
    for (Type& obj : objects)
    {
        obj.updateActors(position);
    }
}


template<class Type>
static void cleanup(PtrList<Type>& objects)
{
    for (Type& obj : objects)
    {
        obj.clear();
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostProcessing::render
(
    vtkMultiProcessController* controller
)
{
    // Some feedback
    if (controller)
    {
        Log << name() << " render (" << controller->GetNumberOfProcesses()
            << " processes)" << endl;
    }
    else
    {
        Log << name() << " render" << endl;
    }

    // Disable any floating point trapping
    // (some low-level rendering functionality does not like it)

    sigFpe::ignore sigFpeHandling; //<- disable in local scope


    // Normal rendering elements
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;

    // Initialise render window
    if (controller || Pstream::master())
    {
        renderer = vtkSmartPointer<vtkRenderer>::New();
        renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

        renderWindow->OffScreenRenderingOn();
        renderWindow->SetSize(output_.width_, output_.height_);

        // Legacy rendering - was deprecated for 8.1.0
        #if (VTK_MAJOR_VERSION < 8) || \
            ((VTK_MAJOR_VERSION == 8) && (VTK_MINOR_VERSION < 2))
        renderWindow->SetAAFrames(10);
        #endif
        renderWindow->SetAlphaBitPlanes(true);
        renderWindow->SetMultiSamples(0);
        // renderWindow->PolygonSmoothingOn();

        renderWindow->AddRenderer(renderer);
    }


    // ---------------------
    #ifdef FOAM_USING_VTK_MPI

    // Multi-process synchronization
    vtkSmartPointer<vtkCompositedSynchronizedRenderers> syncRenderers;
    vtkSmartPointer<vtkSynchronizedRenderWindows> syncWindows;

    if (controller)
    {
        syncRenderers =
            vtkSmartPointer<vtkCompositedSynchronizedRenderers>::New();

        syncWindows =
            vtkSmartPointer<vtkSynchronizedRenderWindows>::New();

        syncWindows->SetRenderWindow(renderWindow);
        syncWindows->SetParallelController(controller);
        syncWindows->SetIdentifier(1);

        // false = Call Render() manually on each process - don't use RMI
        syncWindows->SetParallelRendering(true);

        syncRenderers->SetRenderer(renderer);
        syncRenderers->SetParallelController(controller);
    }
    #endif
    // ---------------------

    scene_.initialise(renderer, output_.name_);

    addGeometryToScene(points_, 0, renderer);
    addGeometryToScene(lines_, 0, renderer);
    addGeometryToScene(surfaces_, 0, renderer);
    addGeometryToScene(text_, 0, renderer);

    while (scene_.loop(renderer))
    {
        const scalar position = scene_.position();

        updateActors(text_, position);
        updateActors(points_, position);
        updateActors(lines_, position);
        updateActors(surfaces_, position);
    }

    // Cleanup
    cleanup(text_);
    cleanup(points_);
    cleanup(lines_);
    cleanup(surfaces_);


    // Instead of relying on the destructor, manually restore the previous
    // SIGFPE state.
    // This is only to avoid compiler complaints about unused variables.

    sigFpeHandling.restore();
}


void Foam::functionObjects::runTimePostProcessing::render
(
    vtkMultiProcessController* controller,
    void* processData
)
{
    reinterpret_cast<runTimePostProcessing*>(processData)->render(controller);
}


void Foam::functionObjects::runTimePostProcessing::render()
{
    #ifdef FOAM_USING_VTK_MPI
    if (parallel_)
    {
        // Create vtkMPIController if MPI is configured,
        // vtkThreadedController otherwise.
        auto ctrl = vtkSmartPointer<vtkMPIController>::New();
        ctrl->Initialize(nullptr, nullptr, 1);

        ctrl->SetSingleMethod(runTimePostProcessing::render, this);
        ctrl->SingleMethodExecute();

        ctrl->Finalize(1);
    }
    else
    #endif
    {
        // Normally we would have a fallback controller like this:

        // if (Pstream::master())
        // {
        //     auto ctrl = vtkSmartPointer<vtkDummyController>::New();
        //     ctrl->Initialize(nullptr, nullptr, 1);
        //
        //     ctrl->SetSingleMethod(runTimePostProcessing::render, this);
        //     ctrl->SingleMethodExecute();
        //
        //     ctrl->Finalize(1);
        // }

        // However, this would prevent us from doing any of our own MPI
        // since this would only be spawned the master.

        // Instead pass in nullptr for the controller and handling
        // logic internally.

        vtkDummyController* dummy = nullptr;
        render(dummy);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostProcessing::runTimePostProcessing
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    output_(),
    parallel_(false),
    scene_(runTime, name),
    points_(),
    lines_(),
    surfaces_(),
    text_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostProcessing::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    #ifdef FOAM_USING_VTK_MPI
    parallel_ = (Pstream::parRun() && dict.getOrDefault("parallel", true));
    #else
    parallel_ = false;
    #endif

    Info<< type() << " " << name() << ": reading post-processing data ("
        << (parallel_ ? "parallel" : "serial") << " rendering)" << endl;

    if (dict.getOrDefault("debug", false))
    {
        debug = 1;
        runTimePostPro::geometryBase::debug = 1;
        runTimePostPro::geometrySurface::debug = 1;
        Info<< "    debugging on" << endl;
    }

    scene_.read(dict);

    const dictionary& outputDict = dict.subDict("output");
    outputDict.readEntry("name", output_.name_);
    outputDict.readEntry("width", output_.width_);
    outputDict.readEntry("height", output_.height_);

    readObjects(dict.subOrEmptyDict("points"), points_);
    readObjects(dict.subOrEmptyDict("lines"), lines_);
    readObjects(dict.subOrEmptyDict("surfaces"), surfaces_);

    const dictionary& textDict = dict.subDict("text");

    for (const entry& dEntry : textDict)
    {
        if (!dEntry.isDict())
        {
            FatalIOErrorInFunction(textDict)
                << textDict.dictName()
                << "text must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const dictionary& objectDict = dEntry.dict();

        text_.append
        (
            new runTimePostPro::text
            (
                *this,
                objectDict,
                scene_.colours()
            )
        );
    }

    return true;
}


bool Foam::functionObjects::runTimePostProcessing::execute()
{
    return true;
}


bool Foam::functionObjects::runTimePostProcessing::write()
{
    render();
    return true;
}


// ************************************************************************* //
