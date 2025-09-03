/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::functionObjects::fieldValues::volFieldValue

Group
    grpFieldFunctionObjects

Description
    Provides a 'volRegion' specialization of the \c fieldValue function object.

    Given a list of user-specified fields and a 'volRegion', a number of
    operations can be performed, such as sums, averages and integrations.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    volFieldValue1
    {
        // Mandatory entries (unmodifiable)
        type            volFieldValue;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields          (<field1> <field2> ... <fieldN>);
        operation       <operationType>;
        regionType      <volRegion>;

        // Optional entries (runtime modifiable)
        postOperation   none;
        weightField     alpha1;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: volFieldValue           | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | Names of operand fields            | wordList | yes | -
      regionType   | Face regionType: see below         | word |  yes  | -
      name         | Name for regionType                | word |  yes  | -
      operation    | Operation type: see below          | word |  yes  | -
      postOperation | Post-operation type: see below    | word |  no   | none
      weightField   | Name of field to apply weighting   | word |  maybe |
      weightFields  | Names of fields to apply weighting | wordList | maybe |
    \endtable

    The inherited entries are elaborated in:
     - \link fieldValue.H \endlink
     - \link volRegion.H \endlink

    Options for the \c regionType entry:
    \plaintable
      cellZone     | requires a 'name' entry to specify the cellZone
      all          | all cells
    \endplaintable

    Options for the \c operation entry:
    \plaintable
      none                 | No operation
      min                  | Minimum
      max                  | Maximum
      sum                  | Sum
      sumMag               | Sum of component magnitudes
      average              | Ensemble average
      volAverage           | Volume weighted average
      volIntegrate         | Volume integral
      CoV                  | Coefficient of variation: standard deviation/mean
      weightedSum          | Weighted sum
      weightedAverage      | Weighted average
      weightedVolAverage   | Weighted volume average
      weightedVolIntegrate | Weighted volume integral
    \endplaintable

    Options for the \c postOperation entry:
    \plaintable
      none          | No additional operation after calculation
      mag           | Component-wise \c mag() after normal operation
      sqrt          | Component-wise \c sqrt() after normal operation
    \endplaintable

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldValues::fieldValue
    - Foam::functionObjects::volRegion

SourceFiles
    volFieldValue.C
    volFieldValueTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volFieldValue_H
#define functionObjects_volFieldValue_H

#include "fieldValue.H"
#include "volRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                        Class volFieldValue Declaration
\*---------------------------------------------------------------------------*/

class volFieldValue
:
    public fieldValue,
    public volRegion
{
public:

    // Public Data Types

        //- Operation type enumeration
        enum operationType
        {
            // Normal operations

            opNone = 0,             //!< No operation
            opMin,                  //!< Minimum value
            opMax,                  //!< Maximum value
            opSum,                  //!< Sum of values
            opSumMag,               //!< Sum of component magnitudes
            opAverage,              //!< Ensemble average
            opVolAverage,           //!< Volume average
            opVolIntegrate,         //!< Volume integral
            opCoV,                  //!< Coefficient of variation

            // Operations returning scalar - bitmask value
            typeScalar = 0x100,

            // Weighted variants
            // Operations using weighting - bitmask value
            typeWeighted = 0x200,

            //! Weighted sum
            opWeightedSum = (opSum | typeWeighted),

            //! Weighted average
            opWeightedAverage = (opAverage | typeWeighted),

            //! Weighted volume average
            opWeightedVolAverage = (opVolAverage | typeWeighted),

            //! Weighted volume integral
            opWeightedVolIntegrate = (opVolIntegrate | typeWeighted),

            // Variants using absolute weighting
            // Operation using mag (eg, for weighting) - bitmask value
            typeAbsolute = 0x400,
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;


        //- Post-operation type enumeration
        enum postOperationType
        {
            postOpNone,   //!< No additional operation after calculation
            postOpMag,    //!< Component-wise mag after normal operation
            postOpSqrt    //!< Component-wise sqrt after normal operation
        };

        //- Operation type names
        static const Enum<postOperationType> postOperationTypeNames_;


protected:

    // Protected Data

        //- Operation to apply to values
        operationType operation_;

        //- Optional post-evaluation operation
        postOperationType postOperation_;

        //- Weight field name(s) - optional
        wordList weightFieldNames_;


    // Protected Member Functions

        //- True if the operation needs the cell-volume
        bool usesVol() const noexcept;

        //- True if the operation variant uses mag
        inline bool is_magOp() const noexcept;

        //- True if the operation variant uses a weight-field
        inline bool is_weightedOp() const noexcept;

        //- True if field is non-empty on any processor.
        inline bool canWeight(const scalarField& fld) const;

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Insert field values into values list
        template<class Type>
        tmp<Field<Type>> getFieldValues
        (
            const word& fieldName,
            const bool mandatory = false
        ) const;

        //- Apply the 'operation' to the values
        template<class Type>
        Type processValues
        (
            const Field<Type>& values,
            const scalarField& V,
            const scalarField& weightField
        ) const;

        //- Helper function to output field values
        label writeAll
        (
            const scalarField& V,
            const scalarField& weightField
        );

        //- Templated helper function to output field values
        template<class Type>
        bool writeValues
        (
            const word& fieldName,
            const scalarField& V,
            const scalarField& weightField
        );

        //- Filter a field according to cellIds
        template<class Type>
        tmp<Field<Type>> filterField(const Field<Type>& field) const;

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("volFieldValue");


    // Constructors

        //- Construct from name, Time and dictionary
        volFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        volFieldValue
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- No copy construct
        volFieldValue(const volFieldValue&) = delete;

        //- No copy assignment
        void operator=(const volFieldValue&) = delete;


    //- Destructor
    virtual ~volFieldValue() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFieldValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volFieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
