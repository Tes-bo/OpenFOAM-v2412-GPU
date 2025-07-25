/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "pointPairs.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Triangulation>
inline Foam::labelPairPair
Foam::pointPairs<Triangulation>::orderPointPair
(
    const labelPair& vA,
    const labelPair& vB
) const
{
    return
    (
        (vA < vB)
      ? labelPairPair(vA, vB)
      : labelPairPair(vB, vA)
    );
}


template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::findPointPair
(
    const labelPair& vA,
    const labelPair& vB
) const
{
    if (vA == vB)
    {
        return false;
    }
    else if (find(orderPointPair(vA, vB)) == end())
    {
        return false;
    }

    return true;
}


template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::insertPointPair
(
    const labelPair& vA,
    const labelPair& vB
)
{
    if (vA == vB)
    {
        return false;
    }

    return insert(orderPointPair(vA, vB));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::pointPairs<Triangulation>::pointPairs(const Triangulation& triangulation)
:
    StorageContainer(),
    triangulation_(triangulation)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::pointPairs<Triangulation>::~pointPairs()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::addPointPair
(
    const labelPair& vA,
    const labelPair& vB
)
{
    return insertPointPair(vA, vB);
}


template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::addPointPair
(
    const labelPair& master,
    const DynamicList<labelPair>& slaves
)
{
    for (const labelPair& slave : slaves)
    {
        addPointPair(master, slave);
    }

    return true;
}


template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::addPointPair
(
    const label vA,
    const label vB
)
{
    const label procNo = Pstream::myProcNo();

    labelPair a(vA, procNo);
    labelPair b(vB, procNo);

    return addPointPair(a, b);
}


template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::isPointPair
(
    const Vertex_handle& vA,
    const Vertex_handle& vB
) const
{
    if (vA->boundaryPoint() && vB->boundaryPoint())
    {
        labelPair a(vA->index(), vA->procIndex());
        labelPair b(vB->index(), vB->procIndex());

        return findPointPair(a, b);
    }

    return false;
}


template<class Triangulation>
inline bool Foam::pointPairs<Triangulation>::isPointPair
(
    const labelPair& vA,
    const labelPair& vB
) const
{
    return findPointPair(vA, vB);
}


template<class Triangulation>
void Foam::pointPairs<Triangulation>::reIndex(const Map<label>& oldToNewIndices)
{
    pointPairs<Triangulation> newTable(triangulation_);

    forAllConstIters(*this, iter)
    {
        labelPairPair e = iter.key();

        labelPair& start = e.first();
        labelPair& end = e.second();

        bool insert = true;

        if (start.second() == Pstream::myProcNo())
        {
            const auto iter2 = oldToNewIndices.cfind(start.first());

            if (iter2.good())
            {
                if (iter2() != -1)
                {
                    start.first() = iter2();
                }
                else
                {
                    insert = false;
                }
            }
        }

        if (end.second() == Pstream::myProcNo())
        {
            const auto iter2 = oldToNewIndices.cfind(end.first());

            if (iter2.good())
            {
                if (iter2() != -1)
                {
                    end.first() = iter2();
                }
                else
                {
                    insert = false;
                }
            }
        }

        if (insert)
        {
            if (e.first() < e.second())
            {
                newTable.insert(e);
            }
            else if (e.first() > e.second())
            {
                newTable.insert(reverse(e));
            }
        }
    }

    this->transfer(newTable);
}


// ************************************************************************* //
