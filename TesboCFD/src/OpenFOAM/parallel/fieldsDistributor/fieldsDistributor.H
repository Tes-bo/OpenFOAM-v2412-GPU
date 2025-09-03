/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    fieldsDistributor

Description
    Common methods/utilities for field decomposers/distributors etc.

SourceFiles
    fieldsDistributorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fieldsDistributor_H
#define Foam_fieldsDistributor_H

#include "bitSet.H"
#include "fileOperation.H"
#include "PtrList.H"
#include "IOobjectList.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fieldsDistributor Declaration
\*---------------------------------------------------------------------------*/

class fieldsDistributor
{
    // Private Methods

    //- Read volume/surface/point/area fields that may or may not exist
    //- on all processors
    template<class BoolListType, class GeoField, class MeshSubsetter>
    static void readFieldsImpl
    (
        refPtr<fileOperation>* readHandlerPtr,  // Can be nullptr
        const BoolListType& haveMeshOnProc,
        const MeshSubsetter* subsetter,         // Can be nullptr
        const typename GeoField::Mesh& mesh,
        IOobjectList& allObjects,
        PtrList<GeoField>& fields,
        const bool deregister
    );


public:

    // Reading helpers

    //- Generic mesh-based field reading
    template<class GeoField>
    static void readField
    (
        const IOobject& io,
        const typename GeoField::Mesh& mesh,
        const label i,
        PtrList<GeoField>& fields
    );

    //- Definition of readField for GeometricFields only
    template<class Type, template<class> class PatchField, class GeoMesh>
    static void readField
    (
        const IOobject& io,
        const typename GeoMesh::Mesh& mesh,
        const label i,
        PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields
    );


    //- Read fields and store on the pointer list
    template
    <
        class Type,
        template<class> class PatchField,
        class GeoMesh
    >
    static void readFields
    (
        const typename GeoMesh::Mesh& mesh,
        const IOobjectList& objects,
        PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
        const bool readOldTime
    );

    //- Read fields and hold on the pointer list
    template<class Mesh, class GeoField>
    static void readFields
    (
        const Mesh& mesh,
        const IOobjectList& objects,
        PtrList<GeoField>& fields
    );


    //- Read volume/surface/point/area fields that may or may not exist
    //- on all processors
    template<class GeoField, class MeshSubsetter>
    static void readFields
    (
        const bitSet& haveMeshOnProc,
        const MeshSubsetter* subsetter,     // Can be nullptr
        const typename GeoField::Mesh& mesh,
        IOobjectList& allObjects,
        PtrList<GeoField>& fields,
        const bool deregister = false
    );

    //- Read volume/surface/point/area fields that may or may not exist
    //- on all processors
    template<class GeoField, class MeshSubsetter>
    static void readFields
    (
        const boolUList& haveMeshOnProc,
        const MeshSubsetter* subsetter,     // Can be nullptr
        const typename GeoField::Mesh& mesh,
        IOobjectList& allObjects,
        PtrList<GeoField>& fields,
        const bool deregister = false
    );

    //- Read volume/surface/point/area fields that may or may not exist
    //- on all processors
    template<class GeoField, class MeshSubsetter>
    static void readFields
    (
        const boolUList& haveMeshOnProc,
        const typename GeoField::Mesh& mesh,
        const autoPtr<MeshSubsetter>& subsetter,
        IOobjectList& allObjects,
        PtrList<GeoField>& fields,
        const bool deregister = false
    );

    //- Read volume/surface/point/area fields that may or may not exist
    //- on all processors
    template<class GeoField, class MeshSubsetter>
    static void readFields
    (
        const boolUList& haveMeshOnProc,
        refPtr<fileOperation>& readHandler,
        const typename GeoField::Mesh& mesh,
        const autoPtr<MeshSubsetter>& subsetter,
        IOobjectList& allObjects,
        PtrList<GeoField>& fields,
        const bool deregister = false
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldsDistributorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
