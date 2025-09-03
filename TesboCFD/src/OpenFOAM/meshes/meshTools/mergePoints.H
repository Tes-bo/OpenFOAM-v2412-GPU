/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Geometric merging of points. See below.

SourceFiles
    mergePoints.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_mergePoints_H
#define Foam_mergePoints_H

#include "scalar.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Detail
{

/*---------------------------------------------------------------------------*\
                   Function Detail::mergePoints Declaration
\*---------------------------------------------------------------------------*/

//- Implementation detail for Foam::mergePoints
//
//  The implementation algorithm is provided with an \em indexer functor
//  to map the subset of points to be merged into the full list.
//  if passed the identityOp this is the same as merging all points.
//
//  \returns Number of point duplicates \em removed by merging.
template<class PointList, class IndexerOp>
label mergePoints
(
    const PointList& points,     //!< The input of all (unmerged) points
    const IndexerOp& indexer,    //!< Indexer for sub-points into all points
    const label nSubPoints,      //!< The number of sub-points
    labelList& pointToUnique,    //!< An old-to-new mapping
    labelList& uniquePoints,     //!< List of unique points from full list
    const scalar mergeTol,       //!< Geometric merge tolerance
    const bool verbose           //!< Debug verbosity
);

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                       Function mergePoints Declarations
\*---------------------------------------------------------------------------*/

//- Calculate merge mapping, preserving the original point order.
//- All points closer/equal mergeTol are to be merged.
//
//  \param[in] points The input (unmerged) points
//  \param[out] pointToUnique  An old-to-new mapping from the original
//      unmerged point index to the index into unique points.
//  \param[out] uniquePoints List of unique points from the original
//      list of unmerged points
//  \param[in] mergeTol Geometric merge tolerance
//  \param[in] verbose Debug verbosity
//
//  \returns Number of point duplicates \em removed by merging.
template<class PointList>
label mergePoints
(
    const PointList& points,
    labelList& pointToUnique,
    labelList& uniquePoints,
    const scalar mergeTol = SMALL,
    const bool verbose = false
);


//- Calculate merge mapping using a subset of points,
//- preserving the original point order.
//- All points closer/equal mergeTol are to be merged.
//
//  \param[in] allPoints The input of all (unmerged) points
//  \param[in] selection The subset of points to consider for merging
//  \param[out] pointToUnique  An old-to-new mapping from the original
//      unmerged point index to the index into unique points.
//  \param[out] uniquePoints List of unique points from the original
//      list of unmerged points
//  \param[in] mergeTol Geometric merge tolerance
//  \param[in] verbose Additional debug verbosity
//
//  \returns Number of point duplicates \em removed by merging.
template<class PointList>
label mergePoints
(
    const PointList& points,
    const labelUList& selection,
    labelList& pointToUnique,
    labelList& uniquePoints,
    const scalar mergeTol = SMALL,
    const bool verbose = false
);


//- Calculate merge mapping, preserving the original point order.
//- All points closer/equal mergeTol are to be merged.
//
//  \param[in] points The input (unmerged) points
//  \param[in] mergeTol Geometric merge tolerance
//  \param[in] verbose Debug verbosity
//  \param[out] pointToUnique An old-to-new mapping from original
//      unmerged point index to the index into unique points.
//
//  \returns Number of \em unique points after merging.
//      Can check against the size of points or pointToUnique to determine
//      if further action (eg, merging, renumbering) is required.
template<class PointList>
label mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique
);


//- Inplace merge points, preserving the original point order.
//- All points closer/equal mergeTol are to be merged.
//
//  \note The list of points must be a concrete List/Field etc.
//  It cannot be an indirect list.
//
//  \param[in,out] points The unmerged points on input, and unique merged
//      points on output
//  \param[in] mergeTol Geometric merge tolerance
//  \param[in] verbose Debug verbosity
//  \param[out] pointToUnique An old-to-new mapping from original
//      unmerged point index to the index into unique points.
//
//  \returns Number of point duplicates \em removed by merging.
//  Can test as true/false to determine further actions.
template<class PointList>
label inplaceMergePoints
(
    PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique
);


//- Inplace merge points a subset of points,
//- preserving the original point order.
//- All points closer/equal mergeTol are to be merged.
//
//  \note The list of points must be a concrete List/Field etc.
//  It cannot be an indirect list.
//
//  \param[in,out] points The unmerged points on input, and unique merged
//      points on output
//  \param[in] selection The subset of points to consider for merging
//  \param[in] mergeTol Geometric merge tolerance
//  \param[in] verbose Debug verbosity
//  \param[out] pointToUnique An old-to-new mapping from original
//      unmerged point index to the index into unique points.
//
//  \returns Number of point duplicates \em removed by merging.
//  Can test as true/false to determine further actions.
template<class PointList>
label inplaceMergePoints
(
    PointList& points,
    const labelUList& selection,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique
);


//- Merge points, preserving the original point order.
//- All points closer/equal mergeTol are to be merged.
//
//  \param[in] points The input (unmerged) points
//  \param[in] mergeTol Geometric merge tolerance
//  \param[in] verbose Debug verbosity
//  \param[out] pointToUnique An old-to-new mapping from original
//      unmerged point index to the index into unique points.
//  \param[out] newPoints The new unique (merged) points.
//
//  \returns Number of point duplicates \em removed by merging.
//  Can test as true/false to determine further actions.
template<class PointList>
label mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointToUnique,
    List<typename PointList::value_type>& newPoints
);



//
// Housekeeping
//

//- Deprecated (MAR-2022) reference point is ignored
//  \deprecated(2022-03) reference point is ignored
template<class PointList>
FOAM_DEPRECATED_FOR(2022-03, "mergePoint() methods without reference point")
inline label mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    typename PointList::value_type origin  /* Used for v2112 and earlier */
)
{
    return Foam::mergePoints(points, mergeTol, verbose, pointMap);
}


//- Deprecated (MAR-2022) reference point is ignored
//  \deprecated(2022-03) reference point is ignored
template<class PointList>
FOAM_DEPRECATED_FOR(2022-03, "mergePoint() methods without reference point")
inline label mergePoints
(
    const PointList& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    List<typename PointList::value_type>& newPoints,
    typename PointList::value_type origin  /* Used for v2112 and earlier */
)
{
    return Foam::mergePoints(points, mergeTol, verbose, pointMap, newPoints);
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mergePoints.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
