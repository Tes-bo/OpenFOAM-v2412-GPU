/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::fieldMinMax

Group
    grpFieldFunctionObjects

Description
    Computes the values and locations of field minima and maxima.
    These are good indicators of calculation performance, e.g. to confirm that
    predicted results are within expected bounds, or how well a case is
    converging.

    Multiple fields can be processed, where for rank > 0 primitives, e.g.
    vectors and tensors, the extrema can be calculated per component, or by
    magnitude.  In addition, spatial location and local processor index are
    included in the output.

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
    fieldMinMax1
    {
        // Mandatory entries (unmodifiable)
        type        fieldMinMax;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        mode        magnitude;
        fields      (<field1> <field2> ... <fieldN>);

        // Optional entries (runtime modifiable)
        location    true;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: fieldMinMax             | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | List of operand fields             | wordList | yes | -
      location     | Write location of the min/max value | bool | no   | true
      mode   | Calculation mode: magnitude or component | word | no | magnitude
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::fieldMinMax

SourceFiles
    fieldMinMax.C
    fieldMinMaxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldMinMax_H
#define functionObjects_fieldMinMax_H

#include "Switch.H"
#include "Enum.H"
#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "vector.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldMinMax Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMax
:
    public fvMeshFunctionObject,
    public writeFile
{

public:

    // Public Enumerations

        enum modeType
        {
            mdMag,      //!< magnitude
            mdCmpt      //!< component
        };


protected:

    // Protected Data

        //- Mode type names
        static const Enum<modeType> modeTypeNames_;

        //- Flag to write location of min/max values
        bool location_;

        //- Mode for min/max - only applicable for ranks > 0
        modeType mode_;

        //- Fields to assess min/max
        volFieldSelection fieldSet_;


    // Protected Member Functions

        //- Helper function to write the output
        template<class Type>
        void output
        (
            const word& fieldName,
            const word& outputName,
            const label minCell,
            const label maxCell,
            const vector& minC,
            const vector& maxC,
            const label minProci,
            const label maxProci,
            const Type& minValue,
            const Type& maxValue
        );


        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Calculate the field min/max for a given field type
        template<class Type>
        void calcMinMaxFieldType
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const word& outputFieldName
        );

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields
        (
            const word& fieldName,
            const modeType& mode
        );


public:

    //- Runtime type information
    TypeName("fieldMinMax");


    // Constructors

        //- Construct from Time and dictionary
        fieldMinMax
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        fieldMinMax(const fieldMinMax&) = delete;

        //- No copy assignment
        void operator=(const fieldMinMax&) = delete;


    //- Destructor
    virtual ~fieldMinMax() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldMinMax
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
