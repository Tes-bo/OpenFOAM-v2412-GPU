/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faPatchMapper.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faPatchMapper::calcAddressing() const
{
    if (directAddrPtr_)
    {
        FatalErrorInFunction
            << "Addressing already calculated"
            << abort(FatalError);
    }

    // Compatibility change  HJ, 12/Aug/2017
    hasUnmapped_ = false;

    directAddrPtr_ = std::make_unique<labelList>(patch_.size(), Foam::zero{});
    auto& addr = *directAddrPtr_;

    // Make a map of old edgeFaces, giving edge index in patch given the new
    // face label next to the patch

    // Create edge index lookup
    Map<label> edgeIndexLookup;

    const labelList& reverseFaceMap = mpm_.reverseFaceMap();

    forAll(oldEdgeFaces_, oefI)
    {
        if (reverseFaceMap[oldEdgeFaces_[oefI]] > -1)
        {
            // Face has survived.  Insert its label under new face index
            edgeIndexLookup.insert(reverseFaceMap[oldEdgeFaces_[oefI]], oefI);
        }
    }

    // Go through new edgeFaces and for each edge try to locate old index
    const labelList& ef = patch_.edgeFaces();

    forAll(ef, efI)
    {
        const auto iter = edgeIndexLookup.cfind(ef[efI]);

        if (iter.good())
        {
            addr[efI] = iter.val();
        }
        else
        {
            // Not found: map from zero
            addr[efI] = 0;

            // Compatibility change  HJ, 12/Aug/2017
            hasUnmapped_ = true;
        }
    }
}


// void Foam::faPatchMapper::clearOut()
// {
//     directAddrPtr_.reset(nullptr);
//     hasUnmapped_ = false;
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faPatchMapper::faPatchMapper
(
    const faPatch& patch,
    const mapPolyMesh& mpm
)
:
    patch_(patch),
    mpm_(mpm),
    sizeBeforeMapping_(patch.size()),
    oldEdgeFaces_(patch.edgeFaces()),
    hasUnmapped_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faPatchMapper::~faPatchMapper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::faPatchMapper::directAddressing() const
{
    if (!directAddrPtr_)
    {
        calcAddressing();
    }

    return *directAddrPtr_;
}


const Foam::labelListList& Foam::faPatchMapper::addressing() const
{
    FatalErrorInFunction
        << "Requested interpolative addressing for a direct mapper."
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::faPatchMapper::weights() const
{
    FatalErrorInFunction
        << "Requested interpolative weights for a direct mapper."
        << abort(FatalError);

    return scalarListList::null();
}


// ************************************************************************* //
