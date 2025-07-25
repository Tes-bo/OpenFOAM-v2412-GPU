/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "Istream.H"
#include "cell.H"
#include "cellModel.H"
#include "IndirectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::cellShape::cellShape() noexcept
:
    labelList(),
    m(nullptr)
{}


inline Foam::cellShape::cellShape
(
    const cellModel& model,
    const labelUList& labels,
    const bool doCollapse
)
:
    labelList(labels),
    m(&model)
{
    if (doCollapse)
    {
        collapse();
    }
}


template<unsigned N>
inline Foam::cellShape::cellShape
(
    const cellModel& model,
    const FixedList<label, N>& labels,
    const bool doCollapse
)
:
    labelList(labels),
    m(&model)
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape
(
    const cellModel& model,
    labelList&& labels,
    const bool doCollapse
)
:
    labelList(std::move(labels)),
    m(&model)
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape
(
    const cellModel::modelType model,
    const labelUList& labels,
    const bool doCollapse
)
:
    labelList(labels),
    m(cellModel::ptr(model))
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape
(
    const cellModel::modelType model,
    labelList&& labels,
    const bool doCollapse
)
:
    labelList(std::move(labels)),
    m(cellModel::ptr(model))
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape
(
    const word& modelName,
    const labelUList& labels,
    const bool doCollapse
)
:
    labelList(labels),
    m(cellModel::ptr(modelName))
{
    if (doCollapse)
    {
        collapse();
    }
}


inline Foam::cellShape::cellShape(Istream& is)
{
    is >> *this;
}


inline Foam::autoPtr<Foam::cellShape> Foam::cellShape::clone() const
{
    return autoPtr<cellShape>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::cellModel& Foam::cellShape::model() const
{
    return *m;
}


inline Foam::label Foam::cellShape::nPoints() const noexcept
{
    return labelList::size();
}


inline Foam::label Foam::cellShape::nEdges() const
{
    return m->nEdges();
}


inline Foam::label Foam::cellShape::nFaces() const
{
    return m->nFaces();
}


inline Foam::pointField Foam::cellShape::points
(
    const UList<point>& meshPoints
) const
{
    return pointField(UIndirectList<point>(meshPoints, *this));
}


inline Foam::labelList Foam::cellShape::meshFaces
(
    const faceList& allFaces,
    const cell& cFaces
) const
{
    // Faces in model order
    faceList localFaces(faces());

    // Do linear match (usually cell shape is low complexity)

    labelList modelToMesh(localFaces.size(), -1);

    forAll(localFaces, i)
    {
        const Foam::face& localF = localFaces[i];

        for (const label meshFacei : cFaces)
        {
            if (allFaces[meshFacei] == localF)
            {
                modelToMesh[i] = meshFacei;
                break;
            }
        }
    }

    return modelToMesh;
}


inline Foam::labelList Foam::cellShape::meshEdges
(
    const edgeList& allEdges,
    const labelList& cEdges
) const
{
    // Edges in model order
    edgeList localEdges(edges());

    // Do linear match (usually cell shape is low complexity)

    labelList modelToMesh(localEdges.size(), -1);

    forAll(localEdges, i)
    {
        const Foam::edge& e = localEdges[i];

        for (const label edgei : cEdges)
        {
            if (allEdges[edgei] == e)
            {
                modelToMesh[i] = edgei;
                break;
            }
        }
    }

    return modelToMesh;
}


inline Foam::face Foam::cellShape::face(const label modelFacei) const
{
    return m->face(modelFacei, *this);
}


inline Foam::faceList Foam::cellShape::faces() const
{
    return m->faces(*this);
}


inline Foam::faceList Foam::cellShape::collapsedFaces() const
{
    const faceList oldFaces(faces());

    faceList newFaces(oldFaces.size());
    label newFacei = 0;

    for (const Foam::face& f : oldFaces)
    {
        Foam::face& newF = newFaces[newFacei];
        newF.resize(f.size());

        label newFp = 0;
        label prevVerti = -1;

        for (const label verti : f)
        {
            if (verti != prevVerti)
            {
                newF[newFp++] = verti;
                prevVerti = verti;
            }
        }

        if ((newFp > 1) && (newF[newFp-1] == newF[0]))
        {
            --newFp;
        }

        if (newFp > 2)
        {
            // Size face and go to next one
            newF.resize(newFp);
            ++newFacei;
        }
    }
    newFaces.resize(newFacei);

    return newFaces;
}


inline Foam::edge Foam::cellShape::edge(const label modelEdgei) const
{
    return m->edge(modelEdgei, *this);
}


inline Foam::edgeList Foam::cellShape::edges() const
{
    return m->edges(*this);
}


inline Foam::point Foam::cellShape::centre(const UList<point>& points) const
{
    return m->centre(*this, points);
}


inline Foam::scalar Foam::cellShape::mag(const UList<point>& points) const
{
    return m->mag(*this, points);
}


inline void Foam::cellShape::reset
(
    const cellModel& model,
    const labelUList& labels,
    const bool doCollapse
)
{
    static_cast<labelList&>(*this) = labels;
    m = &model;

    if (doCollapse)
    {
        collapse();
    }
}


template<unsigned N>
inline void Foam::cellShape::reset
(
    const cellModel& model,
    const FixedList<label, N>& labels,
    const bool doCollapse
)
{
    static_cast<labelList&>(*this) = labels;
    m = &model;

    if (doCollapse)
    {
        collapse();
    }
}


// ************************************************************************* //
