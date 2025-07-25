/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::cuttingSurface

Description
    Constructs a cutting surface through a mesh.

    No attempt at resolving degenerate cases.
    Since the cut faces can be quite ugly, they will often be triangulated.

SourceFiles
    cuttingSurface.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingSurface_H
#define cuttingSurface_H

#include "searchableSurface.H"
#include "cuttingSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cuttingSurface Declaration
\*---------------------------------------------------------------------------*/

class cuttingSurface
:
    public cuttingSurfaceBase
{
    // Private Data

        //- The surface
        const autoPtr<searchableSurface> surfPtr_;


    // Private Member Functions

        //- Determine the cut cells, possibly restricted to a list of cells
        //
        //  \param pointDist [out] For each mesh point, the (signed) distance
        //      to the surface.
        //  \param cellCuts [in,out] On input an empty set (ie, no restriction)
        //      or subsetted cells. On output, the cells cut according to the
        //      planeSides detection.
        void calcCellCuts
        (
            const fvMesh& fvm,
            scalarField& pointDist,
            bitSet& cellCuts
        );


protected:

    // Protected Member Functions

        //- Cut mesh, restricted to a list of cells
        using cuttingSurfaceBase::performCut;

        //- Cut mesh, restricted to a list of cells
        //  Reclaim memory for cellSelectionMask
        virtual void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            bitSet&& cellIdLabels
        );


public:

    // Constructors

        //- Construct from components
        cuttingSurface
        (
            const polyMesh& mesh,
            const word& surfaceType,
            const word& surfaceName
        );

        //- Construct from dictionary
        cuttingSurface
        (
            const word& defaultSurfaceName,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cuttingSurface() = default;


    // Member Functions

        //- The name of the underlying searchableSurface
        const word& surfaceName() const
        {
            return surfPtr_->name();
        }


    // Output

        //- Print information
        void print(Ostream& os, int level=0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
