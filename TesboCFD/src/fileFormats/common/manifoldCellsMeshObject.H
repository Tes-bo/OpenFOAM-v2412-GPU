/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::manifoldCellsMeshObject

Description
    Provides cell-to-faces ('cells()') with duplicate faces removed.

    Useful for postprocessing.

SourceFiles
    manifoldCellsMeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_manifoldCellsMeshObject_H
#define Foam_manifoldCellsMeshObject_H

#include "MeshObject.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class manifoldCellsMeshObject Declaration
\*---------------------------------------------------------------------------*/

class manifoldCellsMeshObject
:
    public MeshObject<polyMesh, UpdateableMeshObject, manifoldCellsMeshObject>
{
    // Private Typedefs

        typedef MeshObject
        <
            polyMesh,
            UpdateableMeshObject,
            manifoldCellsMeshObject
        > MeshObject_type;


    // Private Data

        //- The adjusted cells list
        mutable refPtr<cellList> cellsPtr_;

        //- The number of globally corrected cells (-1 for uninitialized)
        mutable label nCorrected_;


    // Private Member Functions

        //- Calculate cellList with all duplicate faces removed.
        //- Is reference to input mesh if no duplicates
        static refPtr<cellList> filter
        (
            const polyMesh& mesh,
            label& nCellsCorrected
        );

        //- No copy construct
        manifoldCellsMeshObject(const manifoldCellsMeshObject&) = delete;

        //- No copy assignment
        void operator=(const manifoldCellsMeshObject&) = delete;


public:

    //- Declare name of the class, no debug
    ClassNameNoDebug("manifoldCellsMeshObject");


    // Constructors

        //- Construct from mesh
        explicit manifoldCellsMeshObject(const polyMesh& mesh);


    //- Destructor
    virtual ~manifoldCellsMeshObject() = default;


    // Member Functions

        //- Calculate cellList with all duplicate faces removed.
        //- Is reference to input mesh if no duplicates
        static refPtr<cellList> filter(const polyMesh& mesh);

        //- True if any manifold cells detected (globally)
        //- Triggers demand-driven filtering if required.
        bool manifold() const;

        //- Return the (optionally compacted) cell list
        //- Triggers demand-driven filtering if required.
        const cellList& cells() const;

        //- Mesh motion
        virtual bool movePoints()
        {
            return false;
        }

        //- Mesh changes
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
