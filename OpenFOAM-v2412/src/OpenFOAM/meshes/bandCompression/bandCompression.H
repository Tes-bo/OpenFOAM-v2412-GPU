/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    The bandCompression function renumbers the addressing such that the
    band of the matrix is reduced. The algorithm uses a simple search
    through the neighbour list in order of connectivity.

    See http://en.wikipedia.org/wiki/Cuthill-McKee_algorithm

SourceFiles
    bandCompression.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_bandCompression_H
#define Foam_bandCompression_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
template<class T> class CompactListList;

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace meshTools
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Renumber (mesh) addressing to reduce the band of the mesh connectivity,
//- using the Cuthill-McKee algorithm.
//
//  \returns order in which the cells are to be visited (ordered to original)
labelList bandCompression(const polyMesh& mesh);

//- Renumber (mesh) addressing to reduce the band of the matrix,
//- using the Cuthill-McKee algorithm.
//
//  \returns order in which the cells are to be visited (ordered to original)
labelList bandCompression(const CompactListList<label>& addressing);

//- Renumber (mesh) addressing to reduce the band of the matrix,
//- using the Cuthill-McKee algorithm.
//
//  \returns order in which the cells are to be visited (ordered to original)
labelList bandCompression(const labelListList& addressing);

//- Renumber with addressing in losort form (neighbour + start in neighbour),
//- using the Cuthill-McKee algorithm.
//
//  \returns order in which the cells are to be visited (ordered to original)
labelList bandCompression
(
    const labelUList& cellCells,
    const labelUList& offsets
);


} // End namespace meshTools
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Deprecated - prefer meshTools::bandCompression()
//  \deprecated(2022-03) prefer meshTools::bandCompression()
FOAM_DEPRECATED_FOR(2022-03, "meshTools::bandCompression()")
inline labelList bandCompression(const labelListList& cellCellAddressing)
{
    return meshTools::bandCompression(cellCellAddressing);
}

//- Deprecated - prefer meshTools::bandCompression()
//  \deprecated(2022-03) prefer meshTools::bandCompression()
FOAM_DEPRECATED_FOR(2022-03, "meshTools::bandCompression()")
inline labelList bandCompression
(
    const labelUList& cellCells,
    const labelUList& offsets
)
{
    return meshTools::bandCompression(cellCells, offsets);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
