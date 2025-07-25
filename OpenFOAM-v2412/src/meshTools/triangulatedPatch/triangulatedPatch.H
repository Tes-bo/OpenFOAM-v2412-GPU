/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::triangulatedPatch

Description
    Performs a triangulation of a patch to return randomised point locations.

\*---------------------------------------------------------------------------*/

#ifndef Foam_triangulatedPatch_H
#define Foam_triangulatedPatch_H

#include "polyMesh.H"
#include "labelledTri.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class triangulatedPatch Declaration
\*---------------------------------------------------------------------------*/

class triangulatedPatch
{
    // Private Data

        //- Reference to the polyPatch
        const polyPatch& patch_;

        //- Perturbation tolerance to move the point towards the cell centre
        bool perturbTol_;

        //- The polyPatch faces as triangles, the index of each corresponds
        //- to the undecomposed patch face index.
        List<labelledTri> triFace_;

        //- The cumulative triangle area per triangle face,
        //- globally normalised as a 0-1 interval.
        scalarList triWght_;


    // Private Member Functions

        //- Update triangulation
        void update();

        //- Return a random point on the patch based on random number, 0 < c < 1
        bool randomPoint
        (
            Random& rnd,
            const scalar c,
            point& result,
            label& facei,
            label& celli
        ) const;


    // Static Member Functions

        //- Triangulate the patch faces and create addressing
        static void triangulate(const polyPatch& pp, List<labelledTri>& tris);

public:

    //- Constructors

        //- Construct from components
        triangulatedPatch
        (
            const polyPatch& patch,
            const scalar perturbTol
        );

        //- Construct from mesh and patch name
        triangulatedPatch
        (
            const polyMesh& mesh,
            const word& patchName,
            const scalar perturbTol
        );


    //- Destructor
    ~triangulatedPatch() = default;


    // Member Functions

        //- Set a random point on the local patch
        //
        //  \param rnd  random number generator
        //  \param result  the random point
        //  \param facei  index associated with the random point
        //  \param celli  index associated with the random point
        //
        //  \return true  if point has been set
        bool randomLocalPoint
        (
            Random& rnd,
            point& result,
            label& facei,
            label& celli
        ) const;

        //- Set a global random point on the patch
        //
        //  \param  rnd random number generator
        //  \param  result the random point
        //  \param  facei index associated with the random point
        //  \param  celli index associated with the random point
        //
        //  \return true  if point has been set
        bool randomGlobalPoint
        (
            Random& rnd,
            point& result,
            label& facei,
            label& celli
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
