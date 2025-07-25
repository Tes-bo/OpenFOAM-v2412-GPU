/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "faPatchData.H"
#include "dictionary.H"
#include "processorFaPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faPatchData::faPatchData()
:
    ownerPolyPatchId_(-1),
    neighPolyPatchId_(-1),
    ownerProcId_(-1),
    neighProcId_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::faPatchData::name() const noexcept
{
    return name_;
}


Foam::dictionary Foam::faPatchData::dict(const bool withEdgeLabels) const
{
    dictionary patchDict;
    patchDict.add("type", type_);

    if (withEdgeLabels)
    {
        patchDict.add("edgeLabels", edgeLabels_);
    }
    else
    {
        patchDict.add("edgeLabels", labelList());
    }
    patchDict.add("ngbPolyPatchIndex", neighPolyPatchId_);

    if (coupled())
    {
        patchDict.add("myProcNo", ownerProcId_);
        patchDict.add("neighbProcNo", neighProcId_);
    }

    return patchDict;
}


bool Foam::faPatchData::good() const
{
    return (!name_.empty() && !type_.empty());
}


void Foam::faPatchData::clear()
{
    name_.clear();
    type_.clear();

    ownerPolyPatchId_ = -1;
    neighPolyPatchId_ = -1;

    ownerProcId_ = -1;
    neighProcId_ = -1;

    edgeLabels_.clear();
}


void Foam::faPatchData::assign(const faPatch& fap)
{
    clear();

    // Copy information
    name_ = fap.name();
    type_ = fap.type();

    neighPolyPatchId_ = fap.ngbPolyPatchIndex();
    edgeLabels_ = fap.edgeLabels();

    const auto* fapp = isA<processorFaPatch>(fap);
    if (fapp)
    {
        ownerProcId_ = fapp->myProcNo();
        neighProcId_ = fapp->neighbProcNo();
    }
}


bool Foam::faPatchData::assign_coupled(int ownProci, int neiProci)
{
    clear();

    if (ownProci == neiProci)
    {
        return false;
    }

    name_ = processorPolyPatch::newName(ownProci, neiProci);
    type_ = processorFaPatch::typeName;
    ownerProcId_ = ownProci;
    neighProcId_ = neiProci;

    return true;
}


int Foam::faPatchData::matchPatchPair
(
    const labelPair& patchPair
) const noexcept
{
    int ret = 0;
    if (patchPair.first() >= 0 && patchPair.first() == ownerPolyPatchId_)
    {
        ret |= 1;
    }
    if (patchPair.second() >= 0 && patchPair.second() == neighPolyPatchId_)
    {
        ret |= 2;
    }
    return ret;
}


int Foam::faPatchData::comparePatchPair
(
    const labelPair& patchPair
) const noexcept
{
    // As per edge::compare, with validity check

    if (patchPair.first() >= 0 || patchPair.second() >= 0)
    {
        if
        (
            ownerPolyPatchId_ == patchPair.first()
         && neighPolyPatchId_ == patchPair.second()
        )
        {
            return 1;
        }

        if
        (
            ownerPolyPatchId_ == patchPair.second()
         && neighPolyPatchId_ == patchPair.first()
        )
        {
            return -1;
        }
    }

    return 0;
}


// ************************************************************************* //
