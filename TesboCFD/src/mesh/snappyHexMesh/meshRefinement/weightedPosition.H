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

Class
    Foam::weightedPosition

Description
    Wrapper for position + weight to be used in e.g. averaging.

    This avoids the problems when synchronising locations with e.g.
    parallel cyclics. The separation vector only applies to locations
    and not e.g. summed locations. Note that there is no corresponding
    problem for rotational cyclics.

    Typical use might be to e.g. average face centres to points on a patch

    const labelListList& pointFaces = pp.pointFaces();
    const vectorField::subField faceCentres = pp.faceCentres();

    Field<weightedPosition> avgBoundary(pointFaces.size());

    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];
        avgBoundary[pointi].first() = pFaces.size();
        avgBoundary[pointi].second() = sum(pointField(faceCentres, pFaces));
    }
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        weightedPosition::plusEqOp,     // combine op
        pTraits<weightedPosition>::zero,// null value (not used)
        pTraits<weightedPosition>::zero // transform class
    );


SourceFiles
    weightedPosition.C

\*---------------------------------------------------------------------------*/

#ifndef weightedPosition_H
#define weightedPosition_H

#include "Tuple2.H"
#include "point.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class weightedPosition;
class vectorTensorTransform;
class coupledPolyPatch;
class polyMesh;

//- pTraits
template<>
class pTraits<weightedPosition>
{
public:
    typedef weightedPosition cmptType;
    static const weightedPosition zero;
};


/*---------------------------------------------------------------------------*\
                      Class weightedPosition Declaration
\*---------------------------------------------------------------------------*/

class weightedPosition
:
    public Tuple2<scalar, point>
{
public:

    // Constructors

        //- Construct null
        weightedPosition();

        //- Construct from components
        weightedPosition(const scalar s, const point& p);


    // Member Functions

        // Helpers

            //- Get points
            static void getPoints
            (
                const UList<weightedPosition>& in,
                List<point>& out
            );

            //- Set points
            static void setPoints
            (
                const UList<point>& in,
                List<weightedPosition>& out
            );


        //- Summation operator
        static void plusEqOp(weightedPosition& x, const weightedPosition& y);


        // Transformation operators

            void operator()
            (
                const vectorTensorTransform& vt,
                const bool forward,
                List<weightedPosition>& fld
            ) const;

            void operator()
            (
                const vectorTensorTransform& vt,
                const bool forward,
                List<List<weightedPosition>>& flds
            ) const;

            void operator()
            (
                const coupledPolyPatch& cpp,
                Field<weightedPosition>& fld
            ) const;

            template<template<class> class Container>
            void operator()
            (
                const coupledPolyPatch& cpp,
                Container<weightedPosition>& map
            ) const;


        //- Synchronisation for mesh point positions
        static void syncPoints
        (
            const polyMesh& mesh,
            List<weightedPosition>&
        );

        //- Synchronisation for patch point positions
        static void syncPoints
        (
            const polyMesh& mesh,
            const labelUList& meshPoints,
            List<weightedPosition>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "weightedPositionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
