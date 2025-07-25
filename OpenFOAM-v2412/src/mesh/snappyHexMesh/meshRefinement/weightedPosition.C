/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "weightedPosition.H"
#include "vectorTensorTransform.H"
#include "coupledPolyPatch.H"
#include "polyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::weightedPosition Foam::pTraits<Foam::weightedPosition>::zero
(
    scalar(0),
    Foam::point::zero
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::weightedPosition::weightedPosition()
:
    Tuple2<scalar, point>()
{}


Foam::weightedPosition::weightedPosition(const scalar s, const point& p)
:
    Tuple2<scalar, point>(s, p)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::weightedPosition::getPoints
(
    const UList<weightedPosition>& in,
    List<point>& out
)
{
    out.setSize(in.size());
    forAll(in, i)
    {
        out[i] = in[i].second();

        if (mag(in[i].first()) > VSMALL)
        {
            out[i] /= in[i].first();
        }
    }
}


void Foam::weightedPosition::setPoints
(
    const UList<point>& in,
    List<weightedPosition>& out
)
{
    out.setSize(in.size());
    forAll(in, i)
    {
        out[i].second() = out[i].first()*in[i];
    }
}


void Foam::weightedPosition::plusEqOp
(
    weightedPosition& x,
    const weightedPosition& y
)
{
    x.first() += y.first();
    x.second() += y.second();
}


void Foam::weightedPosition::operator()
(
    const vectorTensorTransform& vt,
    const bool forward,
    List<weightedPosition>& fld
) const
{
    pointField pfld;
    getPoints(fld, pfld);

    if (forward)
    {
        pfld = vt.transformPosition(pfld);
    }
    else
    {
        pfld = vt.invTransformPosition(pfld);
    }

    setPoints(pfld, fld);
}


void Foam::weightedPosition::operator()
(
    const vectorTensorTransform& vt,
    const bool forward,
    List<List<weightedPosition>>& flds
) const
{
    for (List<weightedPosition>& fld : flds)
    {
        operator()(vt, forward, fld);
    }
}


void Foam::weightedPosition::operator()
(
    const coupledPolyPatch& cpp,
    Field<weightedPosition>& fld
) const
{
    pointField pfld;
    getPoints(fld, pfld);

    cpp.transformPosition(pfld);

    setPoints(pfld, fld);
}


void Foam::weightedPosition::syncPoints
(
    const polyMesh& mesh,
    List<weightedPosition>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorInFunction << "Size of field " << fld.size()
            << " does not correspond to the number of points in the mesh "
            << mesh.nPoints() << exit(FatalError);
    }

    syncTools::syncPointList
    (
        mesh,
        fld,
        weightedPosition::plusEqOp,     // combine op
        pTraits<weightedPosition>::zero,// null value (not used)
        pTraits<weightedPosition>::zero // transform class
    );
}


void Foam::weightedPosition::syncPoints
(
    const polyMesh& mesh,
    const labelUList& meshPoints,
    List<weightedPosition>& fld
)
{
    if (fld.size() != meshPoints.size())
    {
        FatalErrorInFunction << "Size of field " << fld.size()
            << " does not correspond to the number of points supplied "
            << meshPoints.size() << exit(FatalError);
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        fld,
        weightedPosition::plusEqOp,     // combine op
        pTraits<weightedPosition>::zero,// null value (not used)
        pTraits<weightedPosition>::zero // transform class
    );
}


// ************************************************************************* //
