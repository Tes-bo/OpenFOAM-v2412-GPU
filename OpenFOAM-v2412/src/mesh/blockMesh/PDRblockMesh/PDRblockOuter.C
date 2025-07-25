/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "PDRblock.H"
#include "dictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::PDRblock::outerControl::controlType
>
Foam::PDRblock::outerControl::controlNames_
({
    { controlType::OUTER_NONE,  "none" },
    { controlType::OUTER_EXTEND, "extend" },
    { controlType::OUTER_BOX, "box" },
    { controlType::OUTER_SPHERE, "sphere" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Get a single or a pair of values
template<class T>
static Vector2D<T> getLazyPair(const word& name, const dictionary& dict)
{
    if (token(dict.lookup(name)).isNumber())
    {
        return Vector2D<T>::uniform(dict.get<T>(name));
    }

    return dict.get<Vector2D<T>>(name);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRblock::outerControl::outerControl()
:
    type_(controlType::OUTER_NONE),
    expandType_(expansionType::EXPAND_RATIO),
    onGround_(false),
    relSize_(0,0),
    nCells_(0,0),
    expansion_(1,1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRblock::outerControl::clear()
{
    type_ = controlType::OUTER_NONE;
    expandType_ = expansionType::EXPAND_RATIO;
    onGround_ = false;
    relSize_ = Zero;
    nCells_ = Zero;
    expansion_ = vector2D::uniform(1);
}


void Foam::PDRblock::outerControl::report(Ostream& os) const
{
    if (active())
    {
        os  << "Has outer region: " << controlNames_[type_] << nl
            << " onGround : " << Switch::name(onGround_) << nl
            << "    sizes : " << relSize_ << nl
            << "   nCells : " << nCells_ << nl;
    }
    else
    {
        os  << "No outer region" << nl;
    }
}


bool Foam::PDRblock::outerControl::active() const
{
    return (controlType::OUTER_NONE != type_);
}


bool Foam::PDRblock::outerControl::isSphere() const
{
    return (controlType::OUTER_SPHERE == type_);
}


bool Foam::PDRblock::outerControl::onGround() const
{
    return onGround_;
}


bool Foam::PDRblock::outerControl::onGround(const bool on)
{
    bool old(onGround_);
    onGround_ = on;
    return old;
}


void Foam::PDRblock::outerControl::read(const dictionary& dict)
{
    clear();

    type_ = controlNames_.getOrDefault("type", dict, controlType::OUTER_NONE);
    onGround_ = dict.getOrDefault("onGround", false);

    if (controlType::OUTER_NONE == type_)
    {
        return;
    }

    // Everything else

    nCells_ = getLazyPair<label>("nCells", dict);
    relSize_ = getLazyPair<scalar>("size", dict);

    expandType_ =
        expansionNames_.getOrDefault
        (
            "expansion",
            dict,
            expansionType::EXPAND_RATIO
        );


    if (dict.found("ratios"))
    {
        expansion_ = getLazyPair<scalar>("ratios", dict);
    }
    else
    {
        if (expandType_ != expansionType::EXPAND_UNIFORM)
        {
            expandType_ = expansionType::EXPAND_UNIFORM;
            // Info << "Warning: no 'ratios', use uniform spacing" << nl;
        }
    }

    if (expandType_ == expansionType::EXPAND_UNIFORM)
    {
        expansion_ = vector2D::uniform(1);
    }


    // Errors
    int die = 0;

    if (nCells_.x() <= 1 || nCells_.y() <= 1)
    {
        if (!die++)
        {
            FatalIOErrorInFunction(dict);
        }
        FatalIOError
            << "Too few outer cells: " << nCells_ << nl;
    }

    if (relSize_.x() <= 1 || relSize_.y() <= 1)
    {
        if (!die++)
        {
            FatalIOErrorInFunction(dict);
        }
        FatalIOError
            << "Outer dimensions must be > 1. Had " << relSize_ << nl;
    }

    if (die)
    {
        FatalIOError << nl << exit(FatalIOError);
    }


    // Warnings

    if
    (
        controlType::OUTER_BOX == type_
     || controlType::OUTER_SPHERE == type_
    )
    {
        if (relSize_.x() < 2 || relSize_.y() < 2)
        {
            WarningInFunction
                << "Outer dimensions "
                << relSize_ << " too small for "
                << controlNames_[type_] << " - switching to "
                << controlNames_[controlType::OUTER_EXTEND] << nl;

            type_ = controlType::OUTER_EXTEND;
        }
    }

    if (controlType::OUTER_SPHERE == type_)
    {
        if (relSize_.x() < 3 || relSize_.y() < 3)
        {
            WarningInFunction
                << "Outer dimensions "
                << relSize_ << " too small for "
                << controlNames_[type_] << " - switching to "
                << controlNames_[controlType::OUTER_BOX] << nl;

            type_ = controlType::OUTER_EXTEND;
        }
    }
}


// ************************************************************************* //
