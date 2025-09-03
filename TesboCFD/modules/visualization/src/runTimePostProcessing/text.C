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
#include "text.H"
#include "stringOps.H"
#include "fvMesh.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkCoordinate.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::text::halignType
>
Foam::functionObjects::runTimePostPro::text::halignTypeNames
({
    { halignType::LEFT, "left" },
    { halignType::CENTER, "center" },
    { halignType::CENTER, "centre" },
    { halignType::RIGHT, "right" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::text::text
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    geometryBase(parent, dict, colours),
    string_(dict.get<string>("string")),
    positions_(),
    size_(dict.get<scalar>("size")),
    colour_(nullptr),
    halign_
    (
        halignTypeNames.getOrDefault("halign", dict, halignType::LEFT)
    ),
    bold_(dict.get<bool>("bold")),
    italic_(dict.getOrDefault("italic", false)),
    shadow_(dict.getOrDefault("shadow", false)),
    timeStamp_(dict.getOrDefault("timeStamp", false))
{
    if (!dict.readIfPresent("positions", positions_))
    {
        positions_.resize(1);
        dict.readEntry("position", positions_.front());
    }

    // Additional safety
    if (positions_.empty())
    {
        positions_.emplace_back(0, 0);
    }

    stringOps::inplaceExpand(string_, dict, true, true);

    colour_ = Function1<vector>::NewIfPresent("colour", dict);
    if (!colour_)
    {
        colour_.reset(colours["text"]->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::text::~text()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::text::addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_ || !renderer || !Pstream::master())
    {
        // Add text on master only!
        return;
    }

    DebugInfo << "    Add text: " << string_ << nl;

    // Concatenate string with timeStamp if true
    string str = string_;
    if (timeStamp_)
    {
        str += " " + geometryBase::parent_.mesh().time().timeName();
    }

    const vector textColour = colour_->value(position);

    const scalar textOpacity = opacity(position);

    for (const auto& textPosition : positions_)
    {
        auto actor = vtkSmartPointer<vtkTextActor>::New();

        actor->SetInput(str.c_str());

        vtkTextProperty* prop = actor->GetTextProperty();

        prop->SetFontFamilyToArial();
        prop->SetFontSize(size_);
        prop->SetJustification(int(halign_));
        prop->SetVerticalJustificationToBottom();
        prop->SetBold(bold_);
        prop->SetItalic(italic_);
        prop->SetShadow(shadow_);

        prop->SetColor(textColour[0], textColour[1], textColour[2]);
        prop->SetOpacity(textOpacity);

        // Positioning
        {
            vtkCoordinate* coord = actor->GetPositionCoordinate();

            coord->SetCoordinateSystemToNormalizedViewport();
            coord->SetValue(textPosition.first(), textPosition.second());
        }

        renderer->AddActor2D(actor);
    }
}


void Foam::functionObjects::runTimePostPro::text::updateActors
(
    const scalar position
)
{
    // Do nothing - all handled by addGeometryToScene
}


bool Foam::functionObjects::runTimePostPro::text::clear()
{
    return true;
}


// ************************************************************************* //
