/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "enrichedPatch.H"
#include "primitiveMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcMasterPointFaces() const
{
    if (masterPointFacesPtr_)
    {
        FatalErrorInFunction
            << "Master point face addressing already calculated."
            << abort(FatalError);
    }

    // Note:
    // Master point face addressing lists the master faces for all points
    // in the enriched patch support (if there are no master faces, which is
    // normal, the list will be empty).  The index represents the index of
    // the master face rather than the index from the enriched patch
    // Master face points lists the points of the enriched master face plus
    // points projected into the master face

    Map<DynamicList<label>> mpf(2*meshPoints().size());

    const faceList& ef = enrichedFaces();

    // Add the original face points
    forAll(masterPatch_, facei)
    {
        const face& curFace = ef[facei + slavePatch_.size()];

        for (const label pointi : curFace)
        {
            // Existing or auto-vivify DynamicList
            mpf(pointi).append(facei);
        }
    }

    // Add the projected points which hit the face
    const labelList& slaveMeshPoints = slavePatch_.meshPoints();

    forAll(slavePointFaceHits_, pointi)
    {
        if
        (
            slavePointPointHits_[pointi] < 0
         && slavePointEdgeHits_[pointi] < 0
         && slavePointFaceHits_[pointi].hit()
        )
        {
            // Index of projected point corresponding to this slave point
            const label mergedPointi = pointMergeMap()[slaveMeshPoints[pointi]];

            // Existing or auto-vivify DynamicList
            mpf(mergedPointi).append(slavePointFaceHits_[pointi].hitObject());
        }
    }

    // Re-pack dynamic lists into normal lists

    masterPointFacesPtr_.reset(new Map<labelList>(2*mpf.size()));
    auto& masterPointFaceMap = *masterPointFacesPtr_;

    forAllIters(mpf, mpfIter)
    {
        masterPointFaceMap(mpfIter.key()).transfer(mpfIter.val());
    }
    // Pout<< "masterPointFaceMap: " << masterPointFaceMap << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Map<Foam::labelList>& Foam::enrichedPatch::masterPointFaces() const
{
    if (!masterPointFacesPtr_)
    {
        calcMasterPointFaces();
    }

    return *masterPointFacesPtr_;
}


// ************************************************************************* //
