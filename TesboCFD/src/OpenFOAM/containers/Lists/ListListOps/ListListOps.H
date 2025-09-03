/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

Namespace
    Foam::ListListOps

Description
    Various utility functions to work on Lists of Lists (usually resulting
    from 'gather'ing and combining information from individual processors)

    - combine : \n
        takes (elements of) sublists and appends them into one big list.
    - combineOffset : \n
        similar and also adds offset.

    The access of data is through an AccessOp so that data can be 'gather'ed
    in one go, minimizing communication, and then picked apart and recombined.

    Example:
    \code
        // Assuming myContainer defined which holds all the data I want to
        // transfer (say a pointField and a faceList). myContainer also defines
        // access operators to
        // access the individual elements, say myContainerPoints::operator(),
        // and myContainerFaces::operator()

        List<myContainer> gatheredData(Pstream::nProcs());
        gatheredData[Pstream::myProcNo()] = myContainer(points, faces);

        // Gather data onto master
        Pstream::gatherList(gatheredData);

        // Combine
        pointField combinedPoints
        (
            ListListOps::combine<pointField>
            (
                gatheredData,
                myContainerPoints()
            )
        );

        // Combine and renumber (so combinedFaces indexes combinedPoints)

        // Extract sizes of individual lists
        labelList sizes
        (
            ListListOps::subSizes(gatheredData, myContainerPoints())
        );

        // Renumber using user-defined operator offsetOp<face>()
        faceList combinedFaces
        (
            ListListOps::combineOffset<faceList>
            (
                gatheredData, sizes, myContainerFaces(), offsetOp<face>()
            )
        );
    \endcode

SourceFiles
    ListListOps.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListListOps_H
#define Foam_ListListOps_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Offset operator for ListListOps::combineOffset()
template<class T>
struct offsetOp
{
    T operator()(const T& x, const label offset) const
    {
        return x + offset;
    }
};


/*---------------------------------------------------------------------------*\
                       Namespace ListListOps Declaration
\*---------------------------------------------------------------------------*/

namespace ListListOps
{
    //- Return the sizes of the sub-lists
    template<class T, class AccessOp>
    labelList subSizes(const UList<T>& lists, AccessOp aop = accessOp<T>());

    //- The total size of all sub-lists
    template<class T, class AccessOp>
    label sumSizes(const UList<T>& lists, AccessOp aop = accessOp<T>());

    //- Combines sub-lists into a single list
    template<class AccessType, class T, class AccessOp>
    AccessType combine(const UList<T>& lists, AccessOp aop = accessOp<T>());

    //- Like combine but also offsets sublists based on passed sizes
    template<class AccessType, class T, class AccessOp, class OffsetOp>
    AccessType combineOffset
    (
        const UList<T>& lists,
        const labelUList& offsets,
        AccessOp aop,
        OffsetOp oop = offsetOp<T>()
    );

} // End namespace ListListOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ListListOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
