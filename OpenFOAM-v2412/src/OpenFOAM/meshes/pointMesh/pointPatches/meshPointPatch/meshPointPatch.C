/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "meshPointPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(meshPointPatch, 0);

//- Needs run-time selection table on pointPatch, not facePointPatch
addToRunTimeSelectionTable
(
    pointPatch,
    meshPointPatch,
    dictionary
);

static List<pointConstraint> makeConstraints(const vectorField& normals)
{
    List<pointConstraint> cs(normals.size());

    forAll(cs, i)
    {
        cs[i].applyConstraint(normals[i]);
    }
    return cs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshPointPatch::meshPointPatch
(
    const word& name,
    const labelUList& meshPoints,
    const List<pointConstraint>& constraints,
    const label index,
    const pointBoundaryMesh& bm,
    const word& patchType
)
:
    pointPatch(name, index, bm, word::null, wordList()),
    meshPoints_(meshPoints),
    constraints_(constraints)
{
    if (meshPoints_.size() != constraints_.size())
    {
        FatalErrorInFunction << "patch " << name
            << " size of meshPoints " << meshPoints_.size()
            << " differs from size of constraints " << constraints_.size()
            << exit(FatalError);
    }
}


Foam::meshPointPatch::meshPointPatch
(
    const word& name,
    const labelUList& meshPoints,
    const vectorField& pointNormals,
    const label index,
    const pointBoundaryMesh& bm,
    const word& patchType
)
:
    pointPatch(name, index, bm, word::null, wordList()),
    meshPoints_(meshPoints),
    constraints_(makeConstraints(pointNormals))
{
    if (meshPoints_.size() != pointNormals.size())
    {
        FatalErrorInFunction << "patch " << name
            << " size of meshPoints " << meshPoints_.size()
            << " differs from size of pointNormals " << pointNormals.size()
            << exit(FatalError);
    }
}


Foam::meshPointPatch::meshPointPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointBoundaryMesh& bm,
    const word& patchType
)
:
    pointPatch(name, dict, index, bm),
    meshPoints_(dict.get<labelList>("meshPoints")),
    constraints_(dict.get<List<pointConstraint>>("constraints"))
{}


Foam::meshPointPatch::meshPointPatch
(
    const meshPointPatch& pp,
    const pointBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const labelUList& reversePointMap
)
:
    meshPointPatch
    (
        pp.name(),
        labelList(reversePointMap, labelList(pp.meshPoints(), mapAddressing)),
        List<pointConstraint>(pp.constraints(), mapAddressing),
        index,
        bm,
        pp.type()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshPointPatch::movePoints(PstreamBuffers&, const pointField& p)
{
    localPointsPtr_.reset(nullptr);

    // Recalculate the point normals? Something like
    //if (owner())
    //{
    //    const SubList<face> subFaces
    //    (
    //        mesh.faces(),
    //        mesh.nBoundaryFaces(),
    //        mesh.nInternalFaces()
    //    );
    //    const primitivePatch pp(subFaces, mesh.points());
    //
    //
    //    for (const label pointi : meshPoints())
    //    {
    //        const auto fnd(pp.meshPointMap().find(pointi));
    //        if (fnd)
    //        {
    //            const label patchPointi = fnd();
    //            // Determine point patch equiv
    //
    //        const auto& point
    //
    //

}
    

void Foam::meshPointPatch::updateMesh(PstreamBuffers&)
{
    localPointsPtr_.reset(nullptr);
    pointNormalsPtr_.reset(nullptr);
    // Do what to constraints_? Don't know what the new mesh points are
}


const Foam::pointField& Foam::meshPointPatch::localPoints() const
{
    if (!localPointsPtr_)
    {
        localPointsPtr_.reset
        (
            new pointField
            (
                boundaryMesh().mesh().mesh().points(),
                meshPoints()
            )
        );
    }
    return localPointsPtr_();
}


const Foam::vectorField& Foam::meshPointPatch::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        pointNormalsPtr_.reset(new vectorField(size()));
        vectorField& pointNormals = pointNormalsPtr_();
        forAll(constraints_, i)
        {
            pointNormals[i] = constraints_[i].second();
        }
    }
    return pointNormalsPtr_();
}


void Foam::meshPointPatch::write(Ostream& os) const
{
    pointPatch::write(os);
    meshPoints().writeEntry("meshPoints", os);
    constraints().writeEntry("constraints", os);
}


// ************************************************************************* //
