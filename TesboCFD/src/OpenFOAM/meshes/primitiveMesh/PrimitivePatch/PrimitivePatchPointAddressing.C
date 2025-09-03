/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Description
     Point addressing on the patch: pointEdges and pointFaces.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "DynamicList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcPointEdges() const
{
    DebugInFunction << "Calculating pointEdges" << endl;

    if (pointEdgesPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "pointEdges already calculated"
            << abort(FatalError);
    }

    pointEdgesPtr_.reset(new labelListList(meshPoints().size()));
    auto& pe = *pointEdgesPtr_;

    invertManyToMany(pe.size(), edges(), pe);

    DebugInfo << "    Finished." << endl;
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcPointFaces() const
{
    DebugInFunction
        << "Calculating pointFaces" << endl;

    if (pointFacesPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "pointFaces already calculated"
            << abort(FatalError);
    }

    // Local storage while creating pointFaces
    List<DynamicList<label>> pointFcs(meshPoints().size());

    const List<face_type>& locFcs = localFaces();

    forAll(locFcs, facei)
    {
        for (const label pointi : locFcs[facei])
        {
            pointFcs[pointi].push_back(facei);
        }
    }

    // Copy the lists, recovering content
    pointFacesPtr_.reset(new labelListList(pointFcs.size()));
    auto& pf = *pointFacesPtr_;

    forAll(pointFcs, pointi)
    {
        pf[pointi].transfer(pointFcs[pointi]);
    }

    DebugInfo << "    Finished." << endl;
}


// ************************************************************************* //
