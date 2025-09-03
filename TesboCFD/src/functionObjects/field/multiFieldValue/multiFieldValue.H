/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::functionObjects::fieldValues::multiFieldValue

Group
    grpFieldFunctionObjects

Description
    Computes a selected operation between multiple function objects.

    The operation is applied to all results of each object.

Note
    Each object must generate the same number and type of results.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    multiFieldValue1
    {
        // Mandatory entries
        type        multiFieldValue;
        libs        (fieldFunctionObjects);
        operation   average;

        // List of fieldValue function objects as dictionaries
        functions
        {
            region1
            {
                ...
                // Optional
                resultFields        (field1 field2);
            }
            region2
            {
                ...
                // Optional
                resultFields        (field1 field2);
            }

            ...

            regionN
            {
                ...
                // Optional
                resultFields        (field1 field2);
            }
        }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                         | Type | Reqd | Deflt
      type         | Type name: multiFieldValue          | word |  yes | -
      libs         | Library name: fieldFunctionObjects  | word |  yes | -
      operation    | Operation type to apply to values   | word |  yes | -
      functions    | List of function objects            | dict |  yes | -
    \endtable

    Options for the \c operation entry:
    \plaintable
       sum           | Sum of values
       add           | Add values (same as sum)
       subtract      | Subtract values from first entry
       divide        | Divide first entry by values
       cmptDivide    | Divide first entry by componentwise values
       min           | Minimum value
       max           | Maximum value
       average       | Average value
    \endplaintable

    The \c resultFields entry can be used to set the name of the function object
    result fields to process. If omitted, all available values are employed.

    The inherited entries are elaborated in:
      - \link stateFunctionObject.H \endlink
      - \link writeFile.H \endlink

SourceFiles
    multiFieldValue.C
    multiFieldValueImpl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_multiFieldValue_H
#define Foam_functionObjects_multiFieldValue_H

#include "stateFunctionObject.H"
#include "writeFile.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                       Class multiFieldValue Declaration
\*---------------------------------------------------------------------------*/

class multiFieldValue
:
    public functionObjects::stateFunctionObject,
    public functionObjects::writeFile
{
public:

    // Public Data Types

        //- Operation type enumeration
        enum operationType
        {
            opSum,              //!< Sum of values
            opAdd,              //!< Add values (same as sum)
            opSubtract,         //!< Subtract values from first entry
            opDivide,           //!< Divide first entry by values
            opCmptDivide,       //!< Divide first entry by values componentwise
            opMin,              //!< Minimum value
            opMax,              //!< Maximum value
            opAverage           //!< Average value
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;


private:

    // Private Data

        //- Operation to apply to values
        operationType operation_;

        //- List of function objects
        PtrList<functionObject> functions_;

        //- List of result fields per function object
        List<wordList> resultFields_;


    // Private Member Functions

        //- Templated function to apply the operation.
        //  \return true if Type and resultType are correct
        template<class Type>
        bool applyOperation
        (
            const word& resultType,
            const wordList& names,
            const wordList& entryNames
        );


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader
        (
            const wordList& foNames,
            const List<wordList>& entries,
            const List<wordList>& types,
            Ostream& os
        ) const;


public:

    //- Run-time type information
    TypeName("multiFieldValue");


    // Constructors

        //- Construct from Time and dictionary
        multiFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        multiFieldValue(const multiFieldValue&) = delete;

        //- No copy assignment
        void operator=(const multiFieldValue&) = delete;


    //- Destructor
    virtual ~multiFieldValue() = default;


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
