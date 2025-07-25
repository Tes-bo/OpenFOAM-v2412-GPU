/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::functionObjects::solverInfo

Group
    grpUtilitiesFunctionObjects

Description
    Writes solver information for a list of user-specified fields.

    Information written to file includes:
    - residual fields
    - solver type
    - initial residual
    - final residual
    - number of solver iterations
    - convergence flag

    Operands:
    \table
      Operand      | Type                | Location
      input        | -                   | -
      output file  | dat                 <!--
               --> | postProcessing/\<FO\>/\<time\>/\<file\>(s)
      output field | volScalarField(s)   <!--
               --> | \<time\>/initialResiduals:\<outField\>(s)
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    solverInfo1
    {
        // Mandatory entries (unmodifiable)
        type            solverInfo;
        libs            (utilityFunctionObjects);

        // Mandatory entries (runtime modifiable)
        fields          (<field1> <field2> ... <fieldN>);

        // Optional entries (runtime modifiable)
        writeResidualFields true;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                          | Type | Reqd | Deflt
      type         | Type name: solverInfo                | word | yes  | -
      libs         | Library name: utilityFunctionObjects | word | yes  | -
      fields       | Names of operand fields          | wordList | yes  | -
      writeResidualFields | Flag to write the initial-residual fields <!--
                   -->                                    | bool | no   | false
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

See also
  - Foam::functionObject
  - Foam::functionObjects::fvMeshFunctionObject
  - Foam::functionObjects::writeFile

SourceFiles
    solverInfo.C
    solverInfoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_solverInfo_H
#define functionObjects_solverInfo_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "solverFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class solverInfo Declaration
\*---------------------------------------------------------------------------*/

class solverInfo
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Names of operand fields
        solverFieldSelection fieldSet_;

        //- Names of (result) residual fields
        wordHashSet residualFieldNames_;

        //- Flag to write the initial-residual as a vol field
        bool writeResidualFields_;

        //- Initialisation flag
        bool initialised_;


    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(Ostream& os);

        //- Create and store a residual field on the mesh database
        void createResidualField(const word& fieldName);

        //- Output file header information per primitive type value
        template<class Type>
        void writeFileHeader(Ostream& os, const word& fileName) const;

        //- Initialise a residual field
        template<class Type>
        void initialiseResidualField(const word& fieldName);

        //- Calculate the solver information
        template<class Type>
        void updateSolverInfo(const word& fieldName);


    // Generated Methods

        //- No copy construct
        solverInfo(const solverInfo&) = delete;

        //- No copy assignment
        void operator=(const solverInfo&) = delete;


public:

    //- Runtime type information
    TypeName("solverInfo");


    // Constructors

        //- Construct from Time and dictionary
        solverInfo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~solverInfo() = default;


    // Member Functions

        //- Read solverInfo settings
        virtual bool read(const dictionary&);

        //- Execute solverInfo
        virtual bool execute();

        //- Write solverInfo results
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solverInfoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
