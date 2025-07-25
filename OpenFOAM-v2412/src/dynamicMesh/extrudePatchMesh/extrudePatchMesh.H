/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::extrudePatchMesh

Description
    Mesh at a patch created on the fly. The following entry should be used
    on the field boundary dictionary:

    Example:
    \verbatim
        // New Shell mesh data

        extrudeModel    linearNormal;
        linearNormalCoeffs
        {
            thickness       40e-6;
        }
        nLayers         50;
        expansionRatio  1;
        columnCells      true;

        // Optional feature angle to avoid extruding feature-angles
        featureAngle    45;

        // Patch information
        bottomCoeffs
        {
            name        "bottom";
            type        mappedWall;
            sampleMode  nearestPatchFace;
            samplePatch fixedWalls;
            offsetMode  uniform;
            offset      (0 0 0);
        }

        topCoeffs
        {
            name        "top";
            type        patch;
        }

        sideCoeffs
        {
            name        "side";
            type        empty;
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef Foam_extrudePatchMesh_H
#define Foam_extrudePatchMesh_H

#include "extrudeModel.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class extrudePatchMesh Declaration
\*---------------------------------------------------------------------------*/

class extrudePatchMesh
:
    public fvMesh
{
    // Private Data

        //- Enumeration of patch IDs. Must be 0,1,2,...
        enum patchID
        {
            bottomPatchID = 0,
            topPatchID = 1,
            sidePatchID = 2
        };

        //- Const reference to the patch from which this mesh is extruded
        const polyPatch& extrudedPatch_;

        //- Model dictionary
        dictionary dict_;


    // Private Member Functions

        //- Extrude mesh using polyPatches
        void extrudeMesh(const polyPatchList& regionPatches);

        //- Construct from mesh, patch and dictionary without patches.
        //  Only used internally
        extrudePatchMesh
        (
            const word& regionName,
            const fvMesh& mesh,
            const fvPatch& p,
            const dictionary& dict
        );

public:

    //- Runtime type information
    TypeName("extrudePatchMesh");


    // Constructors

        //- Construct from mesh, patch and dictionary
        extrudePatchMesh
        (
            const fvMesh& mesh,
            const fvPatch& p,
            const dictionary& dict,
            const word& regionName
        );

        //- Construct from mesh, patch, dictionary and new mesh
        //- polyPatch information
        extrudePatchMesh
        (
            const fvMesh& mesh,
            const fvPatch& p,
            const dictionary& dict,
            const word& regionName,
            const polyPatchList& regionPatches
        );

        //- Construct from mesh, patch, dictionary and new mesh
        //- polyPatch information
        extrudePatchMesh
        (
            const fvMesh& mesh,
            const fvPatch& p,
            const dictionary& dict,
            const word& regionName,
            const List<polyPatch*>& regionPatches
        );


    //- Destructor
    virtual ~extrudePatchMesh() = default;


    // Member Functions

        // Access functions

            //- Return region mesh
            const fvMesh& regionMesh() const
            {
                return *this;
            }

            //- Return bottom patch
            const polyPatch& bottomPatch() const
            {
                return this->boundaryMesh()[bottomPatchID];
            }

            //- Return top patch
            const polyPatch& topPatch() const
            {
                 return this->boundaryMesh()[topPatchID];
            }

            //- Return sides patch
            const polyPatch& sidesPatch() const
            {
                 return this->boundaryMesh()[sidePatchID];
            }

            //- Return extruded patch
            const polyPatch& extrudedPatch() const
            {
                return extrudedPatch_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
