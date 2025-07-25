/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::word& Foam::cyclicACMIPolyPatch::nonOverlapPatchName() const
{
    return nonOverlapPatchName_;
}


inline const Foam::polyPatch& Foam::cyclicACMIPolyPatch::nonOverlapPatch() const
{
    // note: use nonOverlapPatchID() as opposed to patch name to initialise
    // demand-driven data

    return this->boundaryMesh()[nonOverlapPatchID()];
}


inline Foam::polyPatch& Foam::cyclicACMIPolyPatch::nonOverlapPatch()
{
    // note: use nonOverlapPatchID() as opposed to patch name to initialise
    // demand-driven data

    return const_cast<polyPatch&>(this->boundaryMesh()[nonOverlapPatchID()]);
}


inline const Foam::scalarField& Foam::cyclicACMIPolyPatch::mask() const
{
    if (owner())
    {
        return srcMask();
    }

    return neighbPatch().tgtMask();
}


inline Foam::scalar Foam::cyclicACMIPolyPatch::tolerance()
{
    return tolerance_;
}


// ************************************************************************* //
