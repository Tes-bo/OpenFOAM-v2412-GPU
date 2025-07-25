/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::processorLODs::cellBox

Description
    Creates the parallel distribution map by describing the source and target
    objects using box shapes.

    A single box is created for the source object, which is then split using
    2x2x2 refinement based on the number of remote target objects that overlap.
    The refinement is local between communicating processor pairs, where the
    refinement continues until a threshold number of remote target objects per
    source box is achieved.

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorLODs_cellBox
#define Foam_processorLODs_cellBox

#include "faceBox.H"
#include "cellList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace processorLODs
{

/*---------------------------------------------------------------------------*\
                           Class cellBox Declaration
\*---------------------------------------------------------------------------*/

class cellBox
:
    public processorLODs::faceBox
{
protected:

    // Protected Data

        //- Reference to the source cell list
        const cellList& srcCells_;

        //- Reference to the target cell list
        const cellList& tgtCells_;


    // Protected Member Functions

        //- The local boundBox associated with given source object
        virtual treeBoundBox calcSrcBox(const label srcObji) const;

        //- The local boundBox associated with given target object
        virtual treeBoundBox calcTgtBox(const label tgtObji) const;



public:

    //- Runtime type information
    TypeName("box");

    // Constructors

        //- Construct from list of points/faces/cells for source and target
        cellBox
        (
            const cellList& srcCells,
            const faceList& srcFaces,
            const UList<point>& srcPoints,
            const cellList& tgtCells,
            const faceList& tgtFaces,
            const UList<point>& tgtPoints,
            const label maxObjectsPerLeaf,
            const label nObjectsOfType,
            const label nRefineIterMax = 100
        );


    //- Destructor
    virtual ~cellBox() = default;


    // Member Functions

        //- Return the parallel distribution map
        //- (usually linear construct order)
        virtual autoPtr<mapDistribute> map
        (
            const mapDistributeBase::layoutTypes constructLayout
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace processorLODs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
