/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "profilingInformation.H"
#include "IOstreams.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::profilingInformation::profilingInformation()
:
    id_(0),
    description_("application::main"),
    parent_(this),
    calls_(0),
    totalTime_(0),
    childTime_(0),
    maxMem_(0),
    active_(false)
{}


Foam::profilingInformation::profilingInformation
(
    profilingInformation *parent,
    const std::string& descr,
    const label id
)
:
    id_(id),
    description_(descr),
    parent_(parent),
    calls_(0),
    totalTime_(0),
    childTime_(0),
    maxMem_(0),
    active_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingInformation::update(const scalar elapsed)
{
    ++calls_;
    totalTime_ += elapsed;

    if (id_ != parent_->id())
    {
        parent_->childTime_ += elapsed;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::profilingInformation::setActive(bool state) const
{
    active_ = state;
}


Foam::Ostream& Foam::profilingInformation::write
(
    Ostream& os,
    const bool offset,
    const scalar elapsedTime,
    const scalar childTimes
) const
{
    os.beginBlock(word("trigger" + Foam::name(id_)));

    os.writeEntry("id",             id_);
    os.writeEntryIfDifferent("parentId", id_, parent().id());
    os.writeEntry("description",    description());
    os.writeEntry("calls",          calls()     + (offset ? 1 : 0));
    os.writeEntry("totalTime",      totalTime() + elapsedTime);
    os.writeEntry("childTime",      childTime() + childTimes);
    os.writeEntryIfDifferent<int>("maxMem", 0, maxMem_);
    os.writeEntry("active",         Switch::name(active()));

    os.endBlock();

    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const profilingInformation& info)
{
    return info.write(os);
}


// ************************************************************************* //
