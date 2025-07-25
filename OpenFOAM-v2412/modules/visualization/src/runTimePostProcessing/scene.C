/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
#include "scene.H"
#include "Constant.H"

// VTK includes
#include "vtkAxesActor.h"
#include "vtkCamera.h"
#include "vtkCaptionActor2D.h"
#include "vtkCubeSource.h"
#include "vtkLightKit.h"
#include "vtkNew.h"
#include "vtkPNGWriter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkWindowToImageFilter.h"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::scene::readCamera
(
    const dictionary& dict
)
{
    nFrameTotal_ = dict.getCheckOrDefault<label>
    (
        "nFrameTotal",
        1,
        labelMinMax::ge(1)
    );

    if
    (
        dict.readCheckIfPresent
        (
            "startPosition",
            startPosition_,
            scalarMinMax::zero_one()
        )
    )
    {
        position_ = startPosition_;
    }

    if (nFrameTotal_ > 1)
    {
        scalar endPosition = dict.getCheckOrDefault<scalar>
        (
            "endPosition",
            1,
            scalarMinMax::zero_one()
        );

        dPosition_ = (endPosition - startPosition_)/scalar(nFrameTotal_ - 1);
    }

    cameraPosition_ = Function1<vector>::New("position", dict);
    cameraFocalPoint_ = Function1<point>::New("focalPoint", dict);
    cameraUp_ = Function1<vector>::New("up", dict);

    dict.readIfPresent("clipBox", clipBox_);
    parallelProjection_ = dict.getOrDefault("parallelProjection", true);
    if (!parallelProjection_)
    {
        cameraViewAngle_ = Function1<scalar>::NewIfPresent("viewAngle", dict);
        if (!cameraViewAngle_)
        {
            cameraViewAngle_.reset
            (
                new Function1Types::Constant<scalar>("viewAngle", 35.0)
            );
        }
    }

    cameraZoom_ = Function1<scalar>::NewIfPresent("zoom", dict);
    if (!cameraZoom_)
    {
        cameraZoom_.reset
        (
            new Function1Types::Constant<scalar>("zoom", 1.0)
        );
    }
}


void Foam::functionObjects::runTimePostPro::scene::readColours
(
    const dictionary& dict
)
{
    const wordList colours(dict.toc());
    for (const word& c : colours)
    {
        colours_.insert(c, Function1<vector>::New(c, dict));
    }
}


void Foam::functionObjects::runTimePostPro::scene::readAxes
(
    const dictionary& dict
)
{
    // Off for backwards compatibility
    axes_.visible = dict.getOrDefault("visible", false);

    if (axes_.visible)
    {
        axes_.colour = dict.getOrDefault<vector>("colour", vector::min);
        axes_.fontColour = dict.getOrDefault<vector>("fontColour", vector::min);
        axes_.fontSize = dict.getOrDefault<scalar>("fontSize", 20);
        axes_.lengthFactor = dict.getOrDefault<scalar>("lengthFactor", 0.1);
    }
}


void Foam::functionObjects::runTimePostPro::scene::setActorVisibility
(
    vtkRenderer* renderer,
    const bool visible
) const
{
    if (!renderer) return;

    vtkActorCollection *actors = renderer->GetActors();
    for (int i = 0; i < actors->GetNumberOfItems(); ++i)
    {
        auto* actor = vtkActor::SafeDownCast(actors->GetItemAsObject(i));
        actor->SetVisibility(visible);
    }
}


void Foam::functionObjects::runTimePostPro::scene::initialise
(
    vtkRenderer* renderer,
    const word& outputName
)
{
    currentFrameI_ = 0;
    position_ = startPosition_;

    outputName_ = outputName;

    if (!renderer) return;


    // Set the background
    const vector bgColour = colours_["background"]->value(position_);

    renderer->SetBackground(bgColour.x(), bgColour.y(), bgColour.z());

    // Apply gradient background if "background2" defined
    if (colours_.contains("background2"))
    {
        const vector bg2Colour = colours_["background2"]->value(position_);

        renderer->GradientBackgroundOn();
        renderer->SetBackground2(bg2Colour.x(), bg2Colour.y(), bg2Colour.z());
    }
    else if (Pstream::parRun())
    {
        // Oddly enough we need a gradient background for parallel rendering,
        // otherwise the colours look quite funny.
        // Doesn't seem to matter if we use SetBackground2() though

        renderer->GradientBackgroundOn();
        renderer->SetBackground2(bgColour.x(), bgColour.y(), bgColour.z());
    }

    // Depth peeling
    renderer->SetUseDepthPeeling(true);
    renderer->SetMaximumNumberOfPeels(4);
    renderer->SetOcclusionRatio(0);

    // Set the camera
    auto camera = vtkSmartPointer<vtkCamera>::New();
    camera->SetParallelProjection(parallelProjection_);
    renderer->SetActiveCamera(camera);

    // Add the lights
    auto lightKit = vtkSmartPointer<vtkLightKit>::New();
    lightKit->AddLightsToRenderer(renderer);

    if (!clipBox_.empty())
    {
        const point& min = clipBox_.min();
        const point& max = clipBox_.max();
        auto clipBox = vtkSmartPointer<vtkCubeSource>::New();
        clipBox->SetXLength(max.x() - min.x());
        clipBox->SetYLength(max.y() - min.y());
        clipBox->SetZLength(max.z() - min.z());
        clipBox->SetCenter
        (
            min.x() + 0.5*(max.x() - min.x()),
            min.y() + 0.5*(max.y() - min.y()),
            min.z() + 0.5*(max.z() - min.z())
        );
        auto clipMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        clipMapper->SetInputConnection(clipBox->GetOutputPort());

        clipBoxActor_ = vtkSmartPointer<vtkActor>::New();
        clipBoxActor_->SetMapper(clipMapper);
        clipBoxActor_->VisibilityOff();
        renderer->AddActor(clipBoxActor_);
    }

    // Add the axes
    if (axes_.visible)
    {
        axes_.actor = vtkSmartPointer<vtkAxesActor>::New();
        renderer->AddActor(axes_.actor);

        axes_.actor->SetTotalLength
        (
            axes_.lengthFactor,
            axes_.lengthFactor,
            axes_.lengthFactor
        );

        auto axisLabel = vtkSmartPointer<vtkTextProperty>::New();
        axisLabel->SetFontFamilyToArial();
        axisLabel->ShadowOff();
        if (axes_.fontColour[0] >= 0)
        {
            const vector& c = axes_.fontColour;
            axisLabel->SetColor(c[0], c[1], c[2]);
        }

        // SetFontSize only works if using SetTextScaleModeToNone - see below
        axisLabel->SetFontSize(axes_.fontSize);

        auto* axisX = axes_.actor->GetXAxisCaptionActor2D();
        auto* axisY = axes_.actor->GetYAxisCaptionActor2D();
        auto* axisZ = axes_.actor->GetZAxisCaptionActor2D();

        axisX->SetCaptionTextProperty(axisLabel);
        axisY->SetCaptionTextProperty(axisLabel);
        axisZ->SetCaptionTextProperty(axisLabel);

        // Note: SetTextScaleModeToNone
        axisX->GetTextActor()->SetTextScaleModeToNone();
        axisY->GetTextActor()->SetTextScaleModeToNone();
        axisZ->GetTextActor()->SetTextScaleModeToNone();

        if (axes_.colour[0] >= 0)
        {
            const vector& c = axes_.colour;
            axes_.actor->GetXAxisShaftProperty()->SetColor(c[0], c[1], c[2]);
            axes_.actor->GetYAxisShaftProperty()->SetColor(c[0], c[1], c[2]);
            axes_.actor->GetZAxisShaftProperty()->SetColor(c[0], c[1], c[2]);
            axes_.actor->GetXAxisTipProperty()->SetColor(c[0], c[1], c[2]);
            axes_.actor->GetYAxisTipProperty()->SetColor(c[0], c[1], c[2]);
            axes_.actor->GetZAxisTipProperty()->SetColor(c[0], c[1], c[2]);
        }
    }
}


void Foam::functionObjects::runTimePostPro::scene::setCamera
(
    vtkRenderer* renderer
) const
{
    if (!renderer) return;

    vtkCamera* camera = renderer->GetActiveCamera();

    if (parallelProjection_)
    {
        // Restore parallel scale to allow application of zoom (later)
        camera->SetParallelScale(1);
    }
    else
    {
        // Restore viewAngle (it might be reset by clipping)
        camera->SetViewAngle(cameraViewAngle_->value(position_));
    }

    const vector up = cameraUp_->value(position_);
    const vector pos = cameraPosition_->value(position_);
    const point focalPoint = cameraFocalPoint_->value(position_);
    const scalar zoom = cameraZoom_->value(position_);

    camera->SetViewUp(up.x(), up.y(), up.z());
    camera->SetPosition(pos.x(), pos.y(), pos.z());
    camera->SetFocalPoint(focalPoint.x(), focalPoint.y(), focalPoint.z());


    // Apply clipping if required
    // Note: possible optimisation - if the camera is static, this only needs
    //       to be done once on initialisation
    if (!clipBox_.empty())
    {
        setActorVisibility(renderer, false);
        clipBoxActor_->VisibilityOn();

        // Call ResetCamera() to fit clip box in view
        renderer->ResetCamera();

        setActorVisibility(renderer, true);
        clipBoxActor_->VisibilityOff();
    }

    // Zoom applied after all other operations
    camera->Zoom(zoom);

    camera->Modified();
}


Foam::string
Foam::functionObjects::runTimePostPro::scene::frameIndexStr() const
{
    string str = Foam::name(currentFrameI_);
    str.insert(0, 4 - str.length(), '0');

    return str;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::scene::scene
(
    const objectRegistry& obr,
    const word& name
)
:
    obr_(obr),
    name_(name),
    colours_(),
    cameraPosition_(nullptr),
    cameraFocalPoint_(nullptr),
    cameraUp_(nullptr),
    cameraViewAngle_(nullptr),
    cameraZoom_(nullptr),
    clipBox_(boundBox::invertedBox),
    clipBoxActor_(),
    axes_(),
    parallelProjection_(true),
    nFrameTotal_(1),
    startPosition_(0),
    position_(0),
    dPosition_(0),
    currentFrameI_(0),
    outputName_("unknown")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::scene::~scene()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::HashPtrTable<Foam::Function1<Foam::vector>, Foam::word>&
Foam::functionObjects::runTimePostPro::scene::colours() const
{
    return colours_;
}


Foam::label Foam::functionObjects::runTimePostPro::scene::frameIndex() const
{
    return currentFrameI_;
}


Foam::scalar Foam::functionObjects::runTimePostPro::scene::position() const
{
    return position_;
}


void Foam::functionObjects::runTimePostPro::scene::read
(
    const dictionary& dict
)
{
    readCamera(dict.subDict("camera"));
    readColours(dict.subDict("colours"));
    readAxes(dict.subOrEmptyDict("axes"));
}


bool Foam::functionObjects::runTimePostPro::scene::loop(vtkRenderer* renderer)
{
    static bool initialised = false;
    setCamera(renderer);

    if (!initialised)
    {
        initialised = true;
        return true;
    }

    if (renderer)
    {

        // Ensure that all objects can be seen without clipping
        // Note: can only be done after all objects have been added!
        renderer->ResetCameraClippingRange();

        // Save image from last iteration
        saveImage(renderer->GetRenderWindow());
    }

    ++currentFrameI_;

    position_ = startPosition_ + currentFrameI_*dPosition_;

    if (currentFrameI_ < nFrameTotal_)
    {
        return true;
    }
    else
    {
        initialised = false;
        return false;
    }
}


void Foam::functionObjects::runTimePostPro::scene::saveImage
(
    vtkRenderWindow* renderWindow
) const
{
    if (!renderWindow)
    {
        return;
    }

    const Time& runTime = obr_.time();

    const fileName fName
    (
        runTime.globalPath()
      / functionObject::outputPrefix
      / name_
      / runTime.timeName()
      / outputName_ + '.' + frameIndexStr() + ".png"
    );

    renderWindow->Render();

    // Set up off-screen rendering
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;

    windowToImageFilter->SetInput(renderWindow);

    //// Add alpha channel for transparency
    // windowToImageFilter->SetInputBufferTypeToRGBA();
    windowToImageFilter->SetInputBufferTypeToRGB();

//    windowToImageFilter->ReadFrontBufferOff();
    windowToImageFilter->Update();

    // Save the image

    if (Pstream::master())
    {
        mkDir(fName.path());

        vtkNew<vtkPNGWriter> writer;
        writer->SetFileName(fName.c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());

        Info<< "    Generating image: " << runTime.relativePath(fName) << endl;
        writer->Write();
    }
}


// ************************************************************************* //
