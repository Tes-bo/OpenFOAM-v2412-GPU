/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::simplifiedMeshes::columnFvMesh

Description
    Generates a 1D column representation of a mesh based on an existing mesh
    and/or fields

SourceFiles
    columnFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedMeshes_columnFvMesh_H
#define simplifiedMeshes_columnFvMesh_H

#include "simplifiedFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedMeshes
{

/*---------------------------------------------------------------------------*\
                        Class columnFvMesh Declaration
\*---------------------------------------------------------------------------*/

class columnFvMeshInfo
{
    // Private Member Functions

        //- Set the patch dictionary entries using the mesh or field files
        bool setPatchEntries(const Time& runTime);

        //- Initialise
        void initialise(const Time& runTime);

        //- Helper function to initialise empty zones
        template<class ZoneMeshType>
        void initialiseZone
        (
            const word& zoneTypeName,
            const fileName& instance,
            ZoneMeshType& zoneType
        );


protected:

    // Protected Data

        //- Region of existing mesh
        const word regionName_;

        //- Location of existing mesh (if present)
        const word localInstance_;

        //- Created from existing mesh (false = from field)
        bool createFromMesh_;

        // Temporary mesh data used during construction of columnFvMesh
        // Note: transferred to columnFvMesh

            //- Points
            pointField points1D_;

            //- Faces
            faceList faces1D_;

            //- Owner addressing
            labelList owner1D_;

            //- Neighbour addressing
            labelList neighbour1D_;


        // Patch data

            //- Dictionary of patch information
            dictionary patchEntries_;

            //- Number of patches with at least 1 local face
            label nPatchWithFace_;


    // Protected Member Functions

        //- Add the patches to the mesh
        void addLocalPatches(fvMesh& mesh) const;

        //- Initialise zones if constructed from mesh
        void initialiseZones(fvMesh& mesh);


public:

    //- Runtime type information
    ClassName("columnFvMeshInfo");

    // Constructor
    columnFvMeshInfo
    (
        const Time& runTime,
        const word& regionName
    );
};


class columnFvMesh
:
    public columnFvMeshInfo,
    public simplifiedFvMesh
{

public:

    //- Runtime type information
    TypeName("columnFvMesh");

    //- Constructor
    columnFvMesh
    (
        const Time& runTime,
        const word& regionName = polyMesh::defaultRegion
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace simplifiedMeshes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "columnFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
