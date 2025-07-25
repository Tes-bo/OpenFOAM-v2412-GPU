/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    Foam::functionObjects::nearWallFields

Group
    grpFieldFunctionObjects

Description
    Samples near-patch volume fields within an input distance range.

    Operands:
    \table
      Operand        | Type             | Location
      input          | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -                | -
      output field   | vol\<Type\>Field | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    nearWallFields1
    {
        // Mandatory entries (unmodifiable)
        type            nearWallFields;
        libs            (fieldFunctionObjects);
        fields
        (
            (<field1> <outField1>)
            (<field2> <outField2>)
        );
        patches         (<patch1> <patch2> ... <patchN>);
        distance        0.01;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd  | Dflt
      type         | Type name: nearWallFields          | word |  yes   | -
      libs         | Library name: fieldFunctionObjects | word |  yes   | -
      fields       | Names of input-output fields | wordHashTable | yes | -
      patches      | Names of patches to sample       | wordList | yes  | -
      distance | Wall-normal distance from patch to sample | scalar | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::nearWallFields

SourceFiles
    nearWallFields.C
    nearWallFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_nearWallFields_H
#define Foam_functionObjects_nearWallFields_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "Tuple2.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class nearWallFields Declaration
\*---------------------------------------------------------------------------*/

class nearWallFields
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        // Read from dictionary

            //- Fields to process (input-name output-name)
            List<Tuple2<word, word>> fieldSet_;

            //- Patches to sample
            labelList patchIDs_;

            //- Distance away from wall
            scalar distance_;

            //- From original field to sampled result
            HashTable<word> fieldMap_;

            //- From resulting back to original field
            HashTable<word> reverseFieldMap_;


        // Calculated addressing

            //- From cell to seed patch faces
            labelListList cellToWalls_;

            //- From cell to tracked end point
            List<List<point>> cellToSamples_;

            //- Map from cell based data back to patch based data
            autoPtr<mapDistribute> getPatchDataMapPtr_;


        // Locally constructed fields

            PtrList<volScalarField> vsf_;
            PtrList<volVectorField> vvf_;
            PtrList<volSphericalTensorField> vSpheretf_;
            PtrList<volSymmTensorField> vSymmtf_;
            PtrList<volTensorField> vtf_;


    // Protected Member Functions

        //- Calculate addressing from cells back to patch faces
        void calcAddressing();

        template<class Type>
        void createFields
        (
            PtrList<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;

        //- Override boundary fields with sampled values
        template<class Type>
        void sampleBoundaryField
        (
            const interpolationCellPoint<Type>& interpolator,
            GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;

        template<class Type>
        void sampleFields
        (
            PtrList<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;


public:

    //- Runtime type information
    TypeName("nearWallFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        nearWallFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        nearWallFields(const nearWallFields&) = delete;

        //- No copy assignment
        void operator=(const nearWallFields&) = delete;


    //- Destructor
    virtual ~nearWallFields() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary& dict);

        //- Calculate the near-wall fields
        virtual bool execute();

        //- Write the near-wall fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nearWallFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
