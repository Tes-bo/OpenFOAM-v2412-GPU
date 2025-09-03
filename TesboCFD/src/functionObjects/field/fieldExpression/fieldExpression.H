/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::functionObjects::fieldExpression

Group
    grpFieldFunctionObjects

Description
    Intermediate class for handling field expression function objects
    (e.g. \c blendingFactor etc.) whereinto a single field is input.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <functionObject>
    {
        // Mandatory and optional entries
        ...

        // Inherited entries
        field             <field>;
        result            <fieldResult>;
        cellZones         ( <cellZone1> <cellZone2> ... <cellZoneN> );
        nLayers           <label>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type  | Reqd | Deflt
      field     | Name of the operand field               | word  | yes  | -
      result    | Name of the output field  | word | no   | \<FO\>(\<field\>)
      cellZones | Names of cellZones                    | wordRes | no   | -
      nLayers   | Number of cell layers around cellZones  | label | no   | 0
    \endtable

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::fieldExpression

SourceFiles
    fieldExpression.C
    fieldExpressionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldExpression_H
#define functionObjects_fieldExpression_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "zoneSubSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class fieldExpression Declaration
\*---------------------------------------------------------------------------*/

class fieldExpression
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Name of field to process
        word fieldName_;

        //- Name of result field
        word resultName_;

        //- Sub-set mesh
        autoPtr<Detail::zoneSubSet> zoneSubSetPtr_;


    // Protected Member Functions

        //- Calculate the components of the field and return true if successful
        virtual bool calc() = 0;

        //- Set the name of result field
        void setResultName(const word& typeName, const word& defaultArg);

        //- Return true if required objects are found
        template<class Type>
        bool foundObject(const word& name, const bool verbose = true) const;


public:

    //- Runtime type information
    TypeName("fieldExpression");


    // Constructors

        //- Construct from Time and dictionary
        fieldExpression
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& fieldName = word::null,
            const word& resultName = word::null
        );

        //- No copy construct
        fieldExpression(const fieldExpression&) = delete;

        //- No copy assignment
        void operator=(const fieldExpression&) = delete;


    //- Destructor
    virtual ~fieldExpression() = default;


    // Member Functions

        //- Read the fieldExpression data
        virtual bool read(const dictionary& dict);

        //- Calculate the result field
        virtual bool execute();

        //- Write the result field
        virtual bool write();

        //- Clear the result field from the objectRegistry
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldExpressionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
