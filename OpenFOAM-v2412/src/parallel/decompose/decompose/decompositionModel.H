/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::decompositionModel

Description
    MeshObject wrapper of decompositionMethod

SourceFiles
    decompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_decompositionModel_H
#define Foam_decompositionModel_H

#include "IOdictionary.H"
#include "MeshObject.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class decompositionModel Declaration
\*---------------------------------------------------------------------------*/

class decompositionModel
:
    public MeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        decompositionModel
    >,
    public IOdictionary
{
    // Private Typedefs

        typedef MeshObject
        <
            polyMesh,
            UpdateableMeshObject,
            decompositionModel
        > MeshObject_type;


    // Private Data

        mutable autoPtr<decompositionMethod> decomposerPtr_;


public:

    // Declare name of the class and its debug switch
    ClassName("decompositionModel");

    //- The canonical name ("decomposeParDict") under which the
    //- MeshObject is registered
    static const word canonicalName;


    // Constructors

        //- Construct from typeName, optional decomposeParDict path/name
        //- or with fallback content
        explicit decompositionModel
        (
            const polyMesh& mesh,
            const fileName& decompDictFile = "",
            const dictionary* fallback = nullptr
        );


    // Selectors

        //- Read and register on mesh,
        //- optionally with alternative decomposeParDict path/name
        //- or with fallback content
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const fileName& decompDictFile = "",
            const dictionary* fallback = nullptr
        );


    // Member Functions

        //- Return demand-driven decomposition method
        decompositionMethod& decomposer() const;


    // UpdateableMeshObject Functions

        virtual bool movePoints()
        {
            return false;
        }

        virtual void updateMesh(const mapPolyMesh&)
        {}


    // Housekeeping

        //- Deprecated(2021-04) compatibility constructor
        //  \deprecated(2021-04)
        FOAM_DEPRECATED_FOR(2021-04, "construct mesh/fileName/dictionary")
        decompositionModel
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const fileName& decompDictFile = ""
        )
        :
            decompositionModel(mesh, decompDictFile, &dict)
        {}

        //- Deprecated(2021-04) compatibility selector
        //  \deprecated(2021-04)
        FOAM_DEPRECATED_FOR(2021-04, "New mesh/fileName/dictionary")
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const fileName& decompDictFile = ""
        )
        {
            return Foam::decompositionModel::New(mesh, decompDictFile, &dict);
        }

        //- Deprecated(2018-08) compatibility method
        //  \deprecated(2018-08) - use IOobject::selectIO directly
        FOAM_DEPRECATED_FOR(2018-08, "IOobject::selectIO")
        static IOobject selectIO
        (
            const IOobject& io,
            const fileName& file,
            const word& name = ""
        )
        {
            return IOobject::selectIO(io, file, name);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
