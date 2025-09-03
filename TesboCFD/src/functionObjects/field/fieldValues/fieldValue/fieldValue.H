/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::functionObjects::fieldValue

Group
    grpFieldFunctionObjects

Description
    Intermediate class for handling field value-based function objects.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <userDefinedSubDictName1>
    {
        // Mandatory and other optional entries
        ...

        // Mandatory (inherited) entries (runtime modifiable)
        fields            (<field1> <field2> ... <fieldN>);

        // Optional (inherited) entries (runtime modifiable)
        writeFields       false;
        scaleFactor       1.0;
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                 | Type     | Req'd | Dflt
      fields      | Names of operand fields     | wordList |  yes  | -
      writeFields | Flag to output field values | bool     |  no   | false
      scaleFactor | Scaling factor              | scalar   |  no   | 1.0
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::fieldValue

SourceFiles
    fieldValue.C
    fieldValueTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldValue_H
#define functionObjects_fieldValue_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldValue Declaration
\*---------------------------------------------------------------------------*/

class fieldValue
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Flag to output field values
        bool writeFields_;

        //- Name of region (patch, zone, etc.)
        word regionName_;

        //- Scaling factor
        scalar scaleFactor_;

        //- Construction dictionary
        dictionary dict_;

        //- Names of operand fields
        wordList fields_;


    // Protected Member Functions

        //- Combine fields from all processor domains into single field
        template<class Type>
        void combineFields(Field<Type>& field);

        //- Combine fields from all processor domains into single field
        template<class Type>
        void combineFields(tmp<Field<Type>>&);


public:

    //- Run-time type information
    TypeName("fieldValue");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fieldValue,
            runTime,
            (
                const word& name,
                const Time& runTime,
                const dictionary& dict
            ),
            (name, runTime, dict)
        );

    // Constructors

        //- Construct from Time and dictionary
        fieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& valueType
        );

        //- Construct from objectRegistry and dictionary
        fieldValue
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const word& valueType
        );

        //- Return a reference to the selected fieldValue
        static autoPtr<fieldValue> New
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool output = true
        );


    //- Destructor
    virtual ~fieldValue() = default;


    // Member Functions

        //- Return the reference to the construction dictionary
        inline const dictionary& dict() const noexcept;

        //- Return the region name
        inline const word& regionName() const noexcept;

        //- Return the list of field names
        inline const wordList& fields() const noexcept;

        //- Return the output field values flag
        inline bool writeFields() const noexcept;

        //- Read from dictionary
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

#include "fieldValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
