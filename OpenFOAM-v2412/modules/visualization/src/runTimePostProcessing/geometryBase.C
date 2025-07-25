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

#include "geometryBase.H"
#include "runTimePostProcessing.H"
#include "Constant.H"

// VTK includes
#include "vtkActor.h"
#include "vtkProperty.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineDebugSwitchWithName(geometryBase, "runTimePostPro::geometryBase", 0);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::geometryBase::renderModeType
>
Foam::functionObjects::runTimePostPro::geometryBase::renderModeTypeNames
({
    { renderModeType::rmFlat, "flat" },
    { renderModeType::rmGouraud, "gouraud" },
    { renderModeType::rmPhong, "phong" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::geometryBase::initialiseActor
(
    vtkActor* actor
) const
{
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetSpecularPower(20);

    switch (renderMode_)
    {
        case rmFlat:
        {
            actor->GetProperty()->SetInterpolationToFlat();
            break;
        }
        case rmGouraud:
        {
            actor->GetProperty()->SetInterpolationToGouraud();
            break;
        }
        case rmPhong:
        {
            actor->GetProperty()->SetInterpolationToPhong();
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryBase::geometryBase
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    parent_(parent),
    name_(dict.dictName()),
    visible_(dict.getOrDefault("visible", true)),
    parallel_
    (
        // User input can only disable parallel here
        #ifdef FOAM_USING_VTK_MPI
        Pstream::parRun() && parent.parallel()
     && dict.getOrDefault("parallel", parent.parallel())
        #else
        false
        #endif
    ),
    renderMode_
    (
        renderModeTypeNames.getOrDefault("renderMode", dict, rmGouraud)
    ),
    opacity_(nullptr),
    colours_(colours)
{
    opacity_ = Function1<scalar>::NewIfPresent("opacity", dict);
    if (!opacity_)
    {
        opacity_.reset(new Function1Types::Constant<scalar>("opacity", 1.0));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryBase::~geometryBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::functionObjects::runTimePostProcessing&
Foam::functionObjects::runTimePostPro::geometryBase::parent() const
{
    return parent_;
}


bool Foam::functionObjects::runTimePostPro::geometryBase::
needsCollective() const
{
    return Pstream::parRun() && (!parent_.parallel() || !parallel_);
}


const Foam::word&
Foam::functionObjects::runTimePostPro::geometryBase::name() const
{
    return name_;
}


Foam::scalar Foam::functionObjects::runTimePostPro::geometryBase::opacity
(
    const scalar position
) const
{
    return opacity_->value(position);
}


const Foam::HashPtrTable<Foam::Function1<Foam::vector>, Foam::word>&
Foam::functionObjects::runTimePostPro::geometryBase::colours() const
{
    return colours_;
}


// ************************************************************************* //
