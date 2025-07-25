/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

const Foam::fvMesh& Foam::VF::raySearchEngine::mesh() const noexcept
{
    return mesh_;
}


const Foam::mapDistribute& Foam::VF::raySearchEngine::map() const
{
    if (!mapPtr_)
    {
        FatalErrorInFunction
            << "mapDistribute has not been set"
            << abort(FatalError);
    }

    return mapPtr_;
}


const Foam::labelList&
Foam::VF::raySearchEngine::compactToGlobal() const noexcept
{
    return compactToGlobal_;
}


const Foam::globalIndex&
Foam::VF::raySearchEngine::globalNumbering() const noexcept
{
    return globalNumbering_;
}


const Foam::labelList& Foam::VF::raySearchEngine::patchIDs() const noexcept
{
    return patchIDs_;
}


const Foam::scalarList& Foam::VF::raySearchEngine::patchAreas() const noexcept
{
    return patchAreas_;
}


Foam::label Foam::VF::raySearchEngine::nParticipatingFaces() const
{
    if (nCoarseFace_ == -1) return nFace_;
    return nCoarseFace_;
}


const Foam::List<Foam::pointField>&
Foam::VF::raySearchEngine::allCf() const noexcept
{
    return allCf_;
}


const Foam::List<Foam::vectorField>&
Foam::VF::raySearchEngine::allSf() const noexcept
{
    return allSf_;
}


const Foam::List<Foam::labelField>&
Foam::VF::raySearchEngine::allAgg() const noexcept
{
    return allAgg_;
}


// ************************************************************************* //
