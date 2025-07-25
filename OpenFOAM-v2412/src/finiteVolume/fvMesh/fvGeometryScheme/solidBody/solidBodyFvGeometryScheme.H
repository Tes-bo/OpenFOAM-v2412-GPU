/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::solidBodyFvGeometryScheme

Description
    Geometry calculation scheme that performs geometry updates only in regions
    where the mesh has changed.

    Example usage in fvSchemes:

    \verbatim
        geometry
        {
            type            solidBody;

            // Optional entries

            // If set to false, update the entire mesh
            partialUpdate   yes;

            // Cache the motion addressing (changed points, faces, cells etc)
            cacheMotion     yes;
        }
    \endverbatim


SourceFiles
    solidBodyFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyFvGeometryScheme_H
#define solidBodyFvGeometryScheme_H

#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class solidBodyFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class solidBodyFvGeometryScheme
:
    public basicFvGeometryScheme
{
    // Private Data

        //- Partial update flag
        bool partialUpdate_;

        //- Cache mesh motion flag
        bool cacheMotion_;

        //- Flag to indicate that the cache has been initialised
        bool cacheInitialised_;

        //- Changed face IDs
        labelList changedFaceIDs_;

        //- Changed patch IDs
        labelList changedPatchIDs_;

        //- Changed cell IDs
        labelList changedCellIDs_;


    // Private Member Functions

        //- Detect what geometry has changed. Return true if anything has.
        bool markChanges
        (
            const pointField& oldPoints,
            const pointField& currPoints,
            bitSet& isChangedPoint,
            bitSet& isChangedFace,
            bitSet& isChangedCell
        ) const;

        //- Set the mesh motion data (point, face IDs)
        void setMeshMotionData();

        //- No copy construct
        solidBodyFvGeometryScheme(const solidBodyFvGeometryScheme&) = delete;

        //- No copy assignment
        void operator=(const solidBodyFvGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("solidBody");


    // Constructors

        //- Construct from mesh
        solidBodyFvGeometryScheme(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~solidBodyFvGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();

        //- Update mesh for topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Calculate geometry quantities using mesh topology and provided
        //- points. If oldPoints provided only does local update. Returns
        //- true if anything changed, false otherwise
        virtual bool updateGeom
        (
            const pointField& points,
            const refPtr<pointField>& oldPoints,    // optional old points
            pointField& faceCentres,
            vectorField& faceAreas,
            pointField& cellCentres,
            scalarField& cellVolumes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
