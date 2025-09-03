/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Patch>
const Foam::scalarListList&
Foam::PrimitivePatchInterpolation<Patch>::faceToPointWeights() const
{
    if (!faceToPointWeightsPtr_)
    {
        makeFaceToPointWeights();
    }

    return *faceToPointWeightsPtr_;
}


template<class Patch>
void Foam::PrimitivePatchInterpolation<Patch>::makeFaceToPointWeights() const
{
    if (faceToPointWeightsPtr_)
    {
        FatalErrorInFunction
            << "Face-to-edge weights already calculated"
            << abort(FatalError);
    }

    const auto& points = patch_.localPoints();
    const auto& faces = patch_.localFaces();

    faceToPointWeightsPtr_.reset(new scalarListList(points.size()));
    auto& weights = *faceToPointWeightsPtr_;

    // get reference to addressing
    const labelListList& pointFaces = patch_.pointFaces();

    forAll(pointFaces, pointi)
    {
        const labelList& curFaces = pointFaces[pointi];

        scalarList& pw = weights[pointi];
        pw.setSize(curFaces.size());

        scalar sumw = 0.0;

        forAll(curFaces, facei)
        {
            pw[facei] =
                1.0/mag(faces[curFaces[facei]].centre(points) - points[pointi]);
            sumw += pw[facei];
        }

        forAll(curFaces, facei)
        {
            pw[facei] /= sumw;
        }
    }
}


template<class Patch>
const Foam::scalarList&
Foam::PrimitivePatchInterpolation<Patch>::faceToEdgeWeights() const
{
    if (!faceToEdgeWeightsPtr_)
    {
        makeFaceToEdgeWeights();
    }

    return *faceToEdgeWeightsPtr_;
}


template<class Patch>
void Foam::PrimitivePatchInterpolation<Patch>::makeFaceToEdgeWeights() const
{
    if (faceToEdgeWeightsPtr_)
    {
        FatalErrorInFunction
            << "Face-to-edge weights already calculated"
            << abort(FatalError);
    }

    const auto& points = patch_.localPoints();
    const auto& faces = patch_.localFaces();
    const edgeList& edges = patch_.edges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    faceToEdgeWeightsPtr_.reset(new scalarList(patch_.nInternalEdges()));
    auto& weights = *faceToEdgeWeightsPtr_;

    forAll(weights, edgei)
    {
        vector P = faces[edgeFaces[edgei][0]].centre(points);
        vector N = faces[edgeFaces[edgei][1]].centre(points);
        vector S = points[edges[edgei].start()];
        vector e = edges[edgei].vec(points);

        scalar alpha =
            -(((N - P)^(S - P))&((N - P)^e))/(((N - P)^e )&((N - P)^e));

        vector E = S + alpha*e;

        weights[edgei] = mag(N - E)/(mag(N - E) + mag(E - P));
    }
}


template<class Patch>
void Foam::PrimitivePatchInterpolation<Patch>::clearWeights()
{
    faceToPointWeightsPtr_.reset(nullptr);
    faceToEdgeWeightsPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Patch>
Foam::PrimitivePatchInterpolation<Patch>::PrimitivePatchInterpolation
(
    const Patch& p
)
:
    patch_(p),
    faceToPointWeightsPtr_(nullptr),
    faceToEdgeWeightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Patch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PrimitivePatchInterpolation<Patch>::faceToPointInterpolate
(
    const Field<Type>& ff
) const
{
    // Check size of the given field
    if (ff.size() != patch_.size())
    {
        FatalErrorInFunction
            << "given field does not correspond to patch. Patch size: "
            << patch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    auto tresult = tmp<Field<Type>>::New(patch_.nPoints(), Zero);
    auto& result = tresult.ref();

    const labelListList& pointFaces = patch_.pointFaces();
    const scalarListList& weights = faceToPointWeights();

    forAll(pointFaces, pointi)
    {
        const labelList& curFaces = pointFaces[pointi];
        const scalarList& w = weights[pointi];

        forAll(curFaces, facei)
        {
            result[pointi] += w[facei]*ff[curFaces[facei]];
        }
    }

    return tresult;
}


template<class Patch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PrimitivePatchInterpolation<Patch>::faceToPointInterpolate
(
    const tmp<Field<Type>>& tff
) const
{
    tmp<Field<Type>> tint = faceToPointInterpolate(tff());
    tff.clear();
    return tint;
}


template<class Patch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PrimitivePatchInterpolation<Patch>::pointToFaceInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != patch_.nPoints())
    {
        FatalErrorInFunction
            << "given field does not correspond to patch. Patch size: "
            << patch_.nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    auto tresult = tmp<Field<Type>>::New(patch_.size(), Zero);
    auto& result = tresult.ref();

    const auto& localFaces = patch_.localFaces();

    forAll(result, facei)
    {
        const labelList& curPoints = localFaces[facei];

        forAll(curPoints, pointi)
        {
            result[facei] += pf[curPoints[pointi]];
        }

        result[facei] /= curPoints.size();
    }

    return tresult;
}


template<class Patch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PrimitivePatchInterpolation<Patch>::pointToFaceInterpolate
(
    const tmp<Field<Type>>& tpf
) const
{
    tmp<Field<Type>> tint = pointToFaceInterpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Patch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PrimitivePatchInterpolation<Patch>::faceToEdgeInterpolate
(
    const Field<Type>& pf
) const
{
    // Check size of the given field
    if (pf.size() != patch_.size())
    {
        FatalErrorInFunction
            << "given field does not correspond to patch. Patch size: "
            << patch_.size() << " field size: " << pf.size()
            << abort(FatalError);
    }

    auto tresult = tmp<Field<Type>>::New(patch_.nEdges(), Zero);
    auto& result = tresult.ref();

    const edgeList& edges = patch_.edges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    const scalarList& weights = faceToEdgeWeights();

    for (label edgei = 0; edgei < patch_.nInternalEdges(); edgei++)
    {
        result[edgei] =
            weights[edgei]*pf[edgeFaces[edgei][0]]
          + (1.0 - weights[edgei])*pf[edgeFaces[edgei][1]];
    }

    for (label edgei = patch_.nInternalEdges(); edgei < edges.size(); edgei++)
    {
        result[edgei] = pf[edgeFaces[edgei][0]];
    }

    return tresult;
}


template<class Patch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PrimitivePatchInterpolation<Patch>::faceToEdgeInterpolate
(
    const tmp<Field<Type>>& tpf
) const
{
    tmp<Field<Type>> tint = faceToEdgeInterpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Patch>
bool Foam::PrimitivePatchInterpolation<Patch>::movePoints()
{
    clearWeights();

    return true;
}


// ************************************************************************* //
