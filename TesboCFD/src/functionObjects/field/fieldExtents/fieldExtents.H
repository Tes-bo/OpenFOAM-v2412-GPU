/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::functionObjects::fieldExtents

Group
    grpFieldFunctionObjects

Description
    Computes the spatial minimum and maximum extents of an input field.

    The extents are derived from the bound box limits after identifying the
    locations where field values exceed the user-supplied threshold value.

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    fieldExtents1
    {
        // Mandatory entries (unmodifiable)
        type                fieldExtents;
        libs                (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields              (<field1> <field2> ... <fieldN>);
        threshold           0.5;

        // Optional entries (runtime modifiable)
        internalField       true;
        referencePosition   (0 0 0);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: fieldExtents            | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      threshold    | Value to identify extents boundary | scalar | yes | -
      fields       | List of operand fields             | wordList | yes | -
      internalField | Flag to process the internal field | bool |  no  | true
      referencePosition | Reference position            | vector | no  | (0 0 0)
      patches  | List of patches to process  | wordList  | no  | \<all patches\>
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the value of \c threshold.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::fieldExtents

SourceFiles
    fieldExtents.C
    fieldExtentsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_fieldExtents_H
#define Foam_functionObjects_fieldExtents_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class fieldExtents Declaration
\*---------------------------------------------------------------------------*/

class fieldExtents
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Flag to write the internal field extents
        bool internalField_;

        //- Threshold value
        scalar threshold_;

        //- Reference position
        point C0_;

        //- Fields to assess
        volFieldSelection fieldSet_;

        //- Patches to assess
        labelList patchIDs_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Return the field mask
        template<class Type>
        tmp<volScalarField> calcMask
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;

        //- Main calculation
        template<class Type>
        void calcFieldExtents
        (
            const word& fieldName,
            const bool calcMag = false
        );


public:

    //- Runtime type information
    TypeName("fieldExtents");


    // Constructors

        //- Construct from Time and dictionary
        fieldExtents
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        fieldExtents(const fieldExtents&) = delete;

        //- No copy assignment
        void operator=(const fieldExtents&) = delete;


    //- Destructor
    virtual ~fieldExtents() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldExtents
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<volScalarField> fieldExtents::calcMask
(
    const GeometricField<scalar, fvPatchField, volMesh>& field
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldExtentsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
