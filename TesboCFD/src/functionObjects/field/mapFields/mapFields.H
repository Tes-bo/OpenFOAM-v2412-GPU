/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::functionObjects::mapFields

Group
    grpFieldFunctionObjects

Description
    Maps input fields from local mesh to secondary mesh at runtime.

    Operands:
    \table
      Operand      | Type                       | Location
      input        | {vol,surface}\<Type\>Field <!--
                                    --> | $FOAM_CASE/\<time\>/\<inpField\>
      output file  | -                          | -
      output field | {vol,surface}\<Type\>Field <!--
                                    --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    mapFields1
    {
        // Mandatory entries (unmodifiable)
        type            mapFields;
        libs            (fieldFunctionObjects);

        // Mandatory (inherited) entries (runtime modifiable)
        fields          (<field1> <field2> ... <fieldN>);
        mapRegion       coarseMesh;
        mapMethod       cellVolumeWeight;
        consistent      true;

        // Optional entries (runtime modifiable)
        // patchMapMethod  direct;  // AMI-related entry

        // When consistent=false

            // Manual: createPatchMap false;
            // patchMap        (<patchSrc> <patchTgt>);
            // cuttingPatches  (<patchTgt1> <patchTgt2> ... <patchTgtN>);

            // Automatic: createPatchMap true;
            // Creates the patchMap and cuttingPatches automatically

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: mapFields               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | Names of operand fields            | wordList |  yes  | -
      mapRegion    | Name of region to map to           | word |  yes  | -
      mapMethod    | Mapping method                     | word |  yes  | -
      consistent   | Mapping meshes have consistent boundaries | bool  | yes | -
      patchMapMethod | Patch mapping method for AMI cases | word |  no  | -
      createPatchMap | Automatically create patchMap and cuttingPatches|  <!--
                 --> | bool | no | -
      patchMap     | Mapping between coincident source and target patches <!--
               --> | wordHashTable | no | -
      cuttingPatches | Target patches cutting the source domain <!--
                 --> | wordList | no | -
    \endtable

    Options for the \c mapMethod entry:
    \verbatim
      direct
      mapNearest
      cellVolumeWeight
      correctedCellVolumeWeight
    \endverbatim

    Options for the \c patchMapMethod entry:
    \verbatim
      nearestFaceAMI
      faceAreaWeightAMI
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::mapFields

SourceFiles
    mapFields.C
    mapFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_mapFields_H
#define functionObjects_mapFields_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshToMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class mapFields Declaration
\*---------------------------------------------------------------------------*/

class mapFields
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Locally cached map region mesh (map to this mesh)
        autoPtr<fvMesh> mapRegionPtr_;

        //- Mesh-to-mesh interpolation
        autoPtr<meshToMesh> interpPtr_;

        //- List of field names to interpolate
        wordRes fieldNames_;


    // Private Member Functions

        //- Helper function to create the mesh-to-mesh interpolation
        void createInterpolation(const dictionary& dict);

        //- Helper function to evaluate constraint patches after mapping
        template<class Type>
        void evaluateConstraintTypes
        (
            GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;

        //- Helper function to map the \<Type\> fields
        template<class Type>
        bool mapFieldType() const;

        //- Helper function to write the \<Type\> fields
        template<class Type>
        bool writeFieldType() const;


public:

    //- Runtime type information
    TypeName("mapFields");


    // Constructors

        //- Construct from Time and dictionary
        mapFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        mapFields(const mapFields&) = delete;

        //- No copy assignment
        void operator=(const mapFields&) = delete;


    //- Destructor
    virtual ~mapFields() = default;


    // Member Functions

        //- Read the mapFields data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the mapFields and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mapFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
