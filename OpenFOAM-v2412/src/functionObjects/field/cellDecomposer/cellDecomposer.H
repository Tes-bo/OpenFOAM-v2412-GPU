/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::functionObjects::cellDecomposer

Group
    grpFieldFunctionObjects

Description
    Maps input fields from local mesh to a secondary mesh at runtime.

    The secondary mesh gets created on-the-fly by decomposing the current mesh.


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
    cellDecomposer1
    {
        // Mandatory entries (unmodifiable)
        type            cellDecomposer;
        libs            (fieldFunctionObjects);

        // Mandatory (inherited) entries (runtime modifiable)
        fields          (<field1> <field2> ... <fieldN>);
        mapRegion       myTetMesh;
        ...
        // Mandatory entries
        // Decompose type: 
        decomposeType   polyhedral;
        // Cell set to decompose
        selectionMode   all;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: cellDecomposer          | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | Names of operand fields            | wordList |  yes  | -
      mapRegion    | Name of region to map to           | word |  yes  | -
      decomposeType | How to decompose cells            | word |  yes  | -
      selectionMode | How to select cells (see fvOption)| word |  yes  | -
    \endtable

    decomposeType:
    - faceCentre : decompose cells into tets using face centre and cell centre.
                   (hex becomes 6*4 tets)
    - faceDiagonal : decompose cells into tets using face diagonal, similar
                     to implicit decomposition inside lagrangian tracking.
                     (hex becomes 6*2 tets)
    - pyramid : keep faces intact but create (polygonal-base) pyramids using
                cell centre (hex becomes 6 pyramids)
    - faceDiagonalQuads : like faceDiagonal but split faces into quads and
                          triangles instead of just triangles
    - polyhedral : like faceDiagonalQuads but only decompose non-hex/prism/tet
                   cells in selected set. Used to convert polyhedral mesh into
                   'simple' mesh.

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

See also
    - Foam::functionObjects::mapFields

SourceFiles
    cellDecomposer.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cellDecomposer_H
#define functionObjects_cellDecomposer_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tetDecomposer;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cellDecomposer Declaration
\*---------------------------------------------------------------------------*/

class cellDecomposer
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Parameter dictionary
        dictionary dict_;

        //- Name of new mesh
        word mapRegion_;

        //- List of field names to interpolate
        wordRes fieldNames_;

        //- Tet decomposer
        autoPtr<tetDecomposer> tetDecompPtr_;

        //- Map from polyMesh to tet-mesh
        autoPtr<mapPolyMesh> mapPtr_;


    // Private Member Functions

        //- Generate mesh
        void makeMesh(const dictionary& dict, const word& name);

        //- Helper function to map the \<Type\> fields
        template<class Type>
        bool mapFieldType() const;

        //- Helper function to write the \<Type\> fields
        template<class Type>
        bool writeFieldType() const;

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const fvMesh& sMesh,
            const labelUList& patchMap,
            const labelUList& cellMap,
            const labelUList& faceMap,
            const bool allowUnmapped
        ) const;


public:


    //- Runtime type information
    TypeName("cellDecomposer");


    // Constructors

        //- Construct from Time and dictionary
        cellDecomposer
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellDecomposer() = default;


    // Member Functions

        //- Read the cellDecomposer data
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellDecomposerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
