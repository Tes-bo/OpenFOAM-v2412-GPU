/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::lduPrimitiveMeshAssembly

Description
    An assembly of lduMatrix that is specific inter-region coupling
    through mapped patches.

SourceFiles
    lduPrimitiveMeshAssembly.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMeshAssembly_H
#define lduPrimitiveMeshAssembly_H

#include "fvMesh.H"
#include "lduPrimitiveMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class lduPrimitiveMeshAssembly Declaration
\*---------------------------------------------------------------------------*/
class lduPrimitiveMeshAssembly
:
    public regIOobject,
    public lduPrimitiveMesh
{
    // Private Data

        //- UPtrList for fvMesh
        UPtrList<lduMesh> meshes_;

        //- Overwrite patchAddr from lduPrimitiveMesh for assembly ldu address
        List<labelList> patchAddr_;

        //- Map from original to asembled (-1 for non-existing)
        labelListList patchMap_;

        //- Map from original to asembled plus extra virtual patch Ids
        labelListList patchLocalToGlobalMap_;

        //- Map internal faces from local to global
        labelListList faceMap_;

        //- Map patch local sub-face to global new internal face
        labelListListList faceBoundMap_;

        //- Map patch local sub-face to nbrCellId
        labelListListList cellBoundMap_;

        //- Map patch local sub-face to local patch face
        labelListListList facePatchFaceMap_;

        //- Offset cells
        labelList cellOffsets_;


    // Private Member Functions

        //- Get size of all meshes
        static label totalSize(const UPtrList<lduMesh>&);

        //- No copy construct
        lduPrimitiveMeshAssembly(const lduPrimitiveMeshAssembly&) = delete;

        //- No copy assignment
        void operator=(const lduPrimitiveMeshAssembly&) = delete;


public:


    // Static Data

        // Declare name of the class and its debug switch
        TypeName("lduPrimitiveMeshAssembly");


    // Constructors

        //- Construct from UPtrList of lduMesh
        lduPrimitiveMeshAssembly
        (
            const IOobject& io,
            const UPtrList<lduMesh>& meshes
        );

        //- Construct from one mesh
        lduPrimitiveMeshAssembly
        (
            const IOobject& io,
            const lduMesh&
        );


    //- Destructor
    virtual ~lduPrimitiveMeshAssembly() = default;


    // Public Member Functions

        //- Update mappings
        template<class Type>
        void update
        (
            UPtrList<GeometricField<Type, fvPatchField, volMesh>>& psis
        );

        //- Size maps
        void updateMaps(const UPtrList<lduMesh>&);

        //- Find nrb mesh Id for mapped patches
        label findNbrMeshId(const polyPatch& pp,const label iMesh) const;


    // Access

        //- Return the object registry
        virtual const objectRegistry& thisDb() const
        {
            return meshes_[0].thisDb();
        }

        //- Return true if thisDb() is a valid DB
        virtual bool hasDb() const
        {
            return true;
        }

        //- Return patch addressing
        virtual const labelUList& patchAddr(const label patchNo) const
        {
            return patchAddr_[patchNo];
        }

        //- Return patchMap
        const labelListList& patchMap() const
        {
            return patchMap_;
        }

        //- Return patchLocalToGlobalMap
        const labelListList& patchLocalToGlobalMap() const
        {
            return patchLocalToGlobalMap_;
        }

        //- Return boundary face map
        const labelListListList& faceBoundMap() const
        {
            return faceBoundMap_;
        }

        //- Return patch local sub-face to nbrCellId map
        const labelListListList& cellBoundMap() const
        {
            return cellBoundMap_;
        }

        //- Return patch local sub-face to local patch face map
        const labelListListList& facePatchFaceMap() const
        {
            return facePatchFaceMap_;
        }

        //- Return cellOffsets
        const labelList& cellOffsets() const
        {
            return cellOffsets_;
        }

        //- Return faceMap
        const labelListList& faceMap() const
        {
            return faceMap_;
        }

        //- Return patchAddr
        const List<labelList>& patchAddr() const
        {
            return patchAddr_;
        }

        //- Return fvMeshes
        const UPtrList<lduMesh>& meshes() const
        {
            return meshes_;
        }

        //- Write Data
        virtual bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduPrimitiveMeshAssemblyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
