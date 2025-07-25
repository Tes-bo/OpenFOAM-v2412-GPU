/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
#include "scalarBar.H"

// #include "doubleVector.H"
// #include "foamVtkTools.H"

// VTK includes
#include "vtkCoordinate.h"
#include "vtkLookupTable.h"
#include "vtkRenderer.h"
#include "vtkScalarBarActor.h"
#include "vtkSmartPointer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::scalarBar::scalarBar()
{
    clear();
}


void Foam::functionObjects::runTimePostPro::scalarBar::clear()
{
    visible_ = true;
    vertical_ = true;
    bold_ = true;
    shadow_ = false;
    italic_ = false;
    titleHack_ = false;
    position_ = {0.8, 0.1};
    size_ = {0.1, 0.5};
    title_ = "";
    fontSize_ = 0;    // 0 == Auto-sizing (defaultFontSize)
    titleSize_ = 0;   // 0 == Auto-sizing (defaultTitleSizeFactor)
    nLabels_ = 5;
    labelFormat_ = "%f";
}


void Foam::functionObjects::runTimePostPro::scalarBar::hide()
{
    visible_ = false;
}


void Foam::functionObjects::runTimePostPro::scalarBar::read
(
    const dictionary& dict
)
{
    clear();

    dict.readIfPresent("visible", visible_);

    if (visible_)
    {
        dict.readIfPresent("vertical", vertical_);
        dict.readIfPresent("bold", bold_);
        dict.readIfPresent("italic", italic_);
        dict.readIfPresent("shadow", shadow_);
        dict.readIfPresent("titleHack", titleHack_);

        if (vertical_)
        {
            size_ = { 0.1, 0.75 };
        }
        else
        {
            size_ = { 0.75, 0.1 };
        }

        dict.readEntry("position", position_);
        dict.readIfPresent("size", size_);
        dict.readEntry("title", title_);
        dict.readIfPresent("fontSize", fontSize_);
        dict.readIfPresent("titleSize", titleSize_);
        dict.readIfPresent("labelFormat", labelFormat_);
        dict.readIfPresent("numberOfLabels", nLabels_);
    }
}


bool Foam::functionObjects::runTimePostPro::scalarBar::add
(
    const vector& textColour,
    vtkRenderer* renderer,
    vtkLookupTable* lut
) const
{
    // Add scalar bar legend

    if (!visible_ || !renderer)
    {
        return false;
    }

    const label fontSizeValue = (fontSize_ ? fontSize_ : defaultFontSize);

    auto sbar = vtkSmartPointer<vtkScalarBarActor>::New();
    sbar->SetLookupTable(lut);
    sbar->SetNumberOfLabels(nLabels_);
    sbar->SetLabelFormat(labelFormat_.c_str());

    /// const vector textColour = colours_["text"]->value(position);

    // Work-around to supply our own scalarbar title
    // - Default scalar bar title text is scales by the scalar bar box
    //   dimensions so if the title is a long string, the text is shrunk to fit
    //   Instead, suppress title and set the title using a vtkTextActor

    vtkSmartPointer<vtkTextActor> titleActor;
    vtkTextProperty* titleProp;

    if (titleHack_)
    {
        // Place the scalar bar title ourselves
        sbar->SetUnconstrainedFontSize(true);

        titleActor = vtkSmartPointer<vtkTextActor>::New();
        titleActor->SetInput(title_.c_str());

        titleProp = titleActor->GetTextProperty();
        titleProp->SetJustificationToCentered();
    }
    else
    {
        // Use the standard scalar bar title
        sbar->SetUnconstrainedFontSize(fontSize_ != 0);
        sbar->SetTitle(title_.c_str());
        titleProp = sbar->GetTitleTextProperty();
    }

    titleProp->SetFontFamilyToArial();

    // Title size was supplied by user (absolute size)
    // or use preset factor (3) of label font size

    if (titleSize_)
    {
        titleProp->SetFontSize(titleSize_);
    }
    else
    {
        // Auto
        titleProp->SetFontSize(defaultTitleSizeFactor*fontSizeValue);

        // Or this??
        // if (!titleHack_) titleProp->SetFontSize(fontSizeValue);
    }

    titleProp->SetJustificationToCentered();
    titleProp->SetVerticalJustificationToBottom();
    titleProp->SetBold(bold_);
    titleProp->SetItalic(italic_);
    titleProp->SetShadow(shadow_);

    titleProp->SetColor(textColour[0], textColour[1], textColour[2]);

    auto labProp = sbar->GetLabelTextProperty();

    labProp->SetColor(textColour[0], textColour[1], textColour[2]);

    if (titleHack_ || fontSize_)
    {
        labProp->SetFontSize(fontSizeValue);
    }
    labProp->ShadowOff();
    labProp->BoldOff();     // or: labProp->SetBold(bold_);
    labProp->ItalicOff();

    // Positioning
    {
        vtkCoordinate* coord = sbar->GetPositionCoordinate();

        coord->SetCoordinateSystemToNormalizedViewport();
        coord->SetValue(position_.first(), position_.second());
    }

    if (vertical_)
    {
        sbar->SetOrientationToVertical();
        sbar->SetTextPositionToSucceedScalarBar();
        sbar->SetWidth(size_.first());
        sbar->SetHeight(size_.second());
        // Standard is sbar->SetBarRatio(0.375);
    }
    else
    {
        sbar->SetOrientationToHorizontal();
        sbar->SetTextPositionToPrecedeScalarBar();

        // Adjustments since not using scalarbar title property
        sbar->SetWidth(size_.first());
        sbar->SetHeight(size_.second());
        // sbar->SetBarRatio(0.5);
        // Standard is sbar->SetBarRatio(0.375);
        // sbar->SetTitleRatio(0.01);
    }

    if (titleActor)
    {
        vtkCoordinate* coord = titleActor->GetPositionCoordinate();

        coord->SetCoordinateSystemToNormalizedViewport();

        coord->SetValue
        (
            position_.first() + (0.5 * sbar->GetWidth()),
            position_.second() + (1.01 * sbar->GetHeight())
        );
    }

    // sbar->DrawFrameOn();
    // sbar->DrawBackgroundOn();
    // sbar->UseOpacityOff();
    // sbar->VisibilityOff();
    sbar->VisibilityOn();

    renderer->AddActor(sbar);

    if (titleActor)
    {
        renderer->AddActor2D(titleActor);
    }

    return true;
}


// ************************************************************************* //
