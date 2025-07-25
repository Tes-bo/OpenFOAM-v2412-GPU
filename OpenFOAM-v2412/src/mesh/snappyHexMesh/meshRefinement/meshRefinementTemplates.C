/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "globalIndex.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::meshRefinement::updateList
(
    const labelList& newToOld,
    const T& nullValue,
    List<T>& elems
)
{
    List<T> newElems(newToOld.size(), nullValue);

    forAll(newElems, i)
    {
        const label oldI = newToOld[i];

        if (oldI >= 0)
        {
            newElems[i] = elems[oldI];
        }
    }

    elems.transfer(newElems);
}


template<class T>
T Foam::meshRefinement::gAverage
(
    const bitSet& isMasterElem,
    const UList<T>& values
)
{
    if (values.size() != isMasterElem.size())
    {
        FatalErrorInFunction
            << "Number of elements in list " << values.size()
            << " does not correspond to number of elements in isMasterElem "
            << isMasterElem.size()
            << exit(FatalError);
    }

    T sum = T(Zero);
    label n = 0;

    forAll(values, i)
    {
        if (isMasterElem.test(i))
        {
            sum += values[i];
            n++;
        }
    }

    reduce(sum, sumOp<T>());
    reduce(n, sumOp<label>());

    if (n > 0)
    {
        return sum/n;
    }
    else
    {
        return pTraits<T>::max;
    }
}


// Compare two lists over all boundary faces
template<class T>
void Foam::meshRefinement::testSyncBoundaryFaceList
(
    const scalar tol,
    const string& msg,
    const UList<T>& faceData,
    const UList<T>& syncedFaceData
) const
{
    const label nBFaces = mesh_.nBoundaryFaces();

    if (faceData.size() != nBFaces || syncedFaceData.size() != nBFaces)
    {
        FatalErrorInFunction
            << "Boundary faces:" << nBFaces
            << " faceData:" << faceData.size()
            << " syncedFaceData:" << syncedFaceData.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        label bFacei = pp.start() - mesh_.nInternalFaces();

        forAll(pp, i)
        {
            const T& data = faceData[bFacei];
            const T& syncData = syncedFaceData[bFacei];

            if (mag(data - syncData) > tol)
            {
                label facei = pp.start()+i;

                FatalErrorInFunction
                    << msg
                    << "patchFace:" << i
                    << " face:" << facei
                    << " fc:" << mesh_.faceCentres()[facei]
                    << " patch:" << pp.name()
                    << " faceData:" << data
                    << " syncedFaceData:" << syncData
                    << " diff:" << mag(data - syncData)
                    << abort(FatalError);
            }

            bFacei++;
        }
    }
}


// Print list sorted by coordinates. Used for comparing non-parallel v.s.
// parallel operation
template<class T>
void Foam::meshRefinement::collectAndPrint
(
    const UList<point>& points,
    const UList<T>& data
)
{
    const globalIndex globalPoints(points.size());

    pointField allPoints(globalPoints.gather(points));
    List<T> allData(globalPoints.gather(data));


    scalarField magAllPoints(mag(allPoints-point(-0.317, 0.117, 0.501)));

    labelList visitOrder(sortedOrder(magAllPoints));
    for (const label allPointi : visitOrder)
    {
        Info<< allPoints[allPointi] << " : " << allData[allPointi]
            << endl;
    }
}


template<class GeoField>
void Foam::meshRefinement::addPatchFields
(
    fvMesh& mesh,
    const word& patchFieldType
)
{
    // Field order is irrelevant...
    for (GeoField& fld : mesh.objectRegistry::objects<GeoField>())
    {
        auto& bfld = fld.boundaryFieldRef();

        const label newPatchi = bfld.size();
        bfld.resize(newPatchi+1);

        bfld.set
        (
            newPatchi,
            GeoField::Patch::New
            (
                patchFieldType,
                mesh.boundary()[newPatchi],
                fld.internalField()
            )
        );
    }
}


template<class GeoField>
void Foam::meshRefinement::reorderPatchFields
(
    fvMesh& mesh,
    const labelList& oldToNew
)
{
    // Field order is irrelevant...
    for (GeoField& fld : mesh.objectRegistry::objects<GeoField>())
    {
        fld.boundaryFieldRef().reorder(oldToNew);
    }
}


template<class EnumContainer>
int Foam::meshRefinement::readFlags
(
    const EnumContainer& namedEnum,
    const wordList& words
)
{
    int flags = 0;

    for (const word& w : words)
    {
        flags |= namedEnum[w];
    }

    return flags;
}


template<class Type>
void Foam::meshRefinement::weightedSum
(
    const polyMesh& mesh,
    const bitSet& isMasterEdge,
    const labelList& meshPoints,
    const edgeList& edges,
    const scalarField& edgeWeights,
    const Field<Type>& pointData,
    Field<Type>& sum
)
{
    if
    (
        edges.size() != isMasterEdge.size()
     || edges.size() != edgeWeights.size()
     || meshPoints.size() != pointData.size()
    )
    {
        FatalErrorInFunction
            << "Inconsistent sizes for edge or point data:"
            << " isMasterEdge:" << isMasterEdge.size()
            << " edgeWeights:" << edgeWeights.size()
            << " edges:" << edges.size()
            << " pointData:" << pointData.size()
            << " meshPoints:" << meshPoints.size()
            << abort(FatalError);
    }

    sum.setSize(meshPoints.size());
    sum = Type(Zero);

    forAll(edges, edgeI)
    {
        if (isMasterEdge.test(edgeI))
        {
            const edge& e = edges[edgeI];

            scalar eWeight = edgeWeights[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            sum[v0] += eWeight*pointData[v1];
            sum[v1] += eWeight*pointData[v0];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        sum,
        plusEqOp<Type>(),
        Type(Zero)          // null value
    );
}


template<class Type>
Type Foam::meshRefinement::get
(
    const dictionary& dict,
    const word& keyword,
    const bool noExit,
    enum keyType::option matchOpt,
    const Type& deflt
)
{
    // If noExit, emit FatalIOError message without exiting
    IOobjectOption::readOption readOpt
    (
        noExit
      ? IOobjectOption::READ_IF_PRESENT : IOobjectOption::MUST_READ
    );

    Type val(deflt);

    if (!dict.readEntry(keyword, val, matchOpt, readOpt))
    {
        FatalIOErrorInFunction(dict)
            << "Entry '" << keyword << "' not found in dictionary "
            << dict.name() << nl;
    }

    return val;
}


// ************************************************************************* //
