/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::functionObjects::processorField

Group
    grpFieldFunctionObjects

Description
    Writes a scalar field whose value is the local processor ID.  The output
    field name is \c processorID.

    Operands:
    \table
      Operand        | Type           | Location
      input          | volScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    processorField1
    {
        // Mandatory entries (unmodifiable)
        type        processorField;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: processorField          | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func processorField
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::processorField

SourceFiles
    processorField.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_processorField_H
#define functionObjects_processorField_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class processorField Declaration
\*---------------------------------------------------------------------------*/

class processorField
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("processorField");


    // Constructors

        //- Construct from Time and dictionary
        processorField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        processorField(const processorField&) = delete;

        //- No copy assignment
        void operator=(const processorField&) = delete;


    //- Destructor
    virtual ~processorField() = default;


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary&);

        //- Calculate the processorID field
        virtual bool execute();

        //- Write the processorID field
        virtual bool write();

        //- Update for changes of mesh
        //  The base implementation is a no-op.
        virtual void updateMesh(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
