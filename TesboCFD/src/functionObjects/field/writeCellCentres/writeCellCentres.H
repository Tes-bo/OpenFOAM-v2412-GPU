/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::writeCellCentres

Group
    grpFieldFunctionObjects

Description
    Writes the cell-centres volVectorField and the three component fields as
    \c volScalarFields.

    Operands:
    \table
      Operand          | Type           | Location
      input            | -              | -
      output file      | -              | -
      output field 1   | volVectorField | $FOAM_CASE/\<time\>/C
      output field 2   | volScalarField | $FOAM_CASE/\<time\>/{Cx,Cy,Cz}
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    writeCellCentres
    {
        // Mandatory entries (unmodifiable)
        type        writeCellCentres;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                        | Type | Req'd | Dflt
      type       | Type name: writeCellCentres        | word |  yes  | -
      libs       | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func writeCellCentres
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::writeCellCentres

SourceFiles
    writeCellCentres.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeCellCentres_H
#define functionObjects_writeCellCentres_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class writeCellCentres Declaration
\*---------------------------------------------------------------------------*/

class writeCellCentres
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("writeCellCentres");


    // Constructors

        //- Construct from Time and dictionary
        writeCellCentres
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        writeCellCentres(const writeCellCentres&) = delete;

        //- No copy assignment
        void operator=(const writeCellCentres&) = delete;


    //- Destructor
    virtual ~writeCellCentres() = default;


    // Member Functions

        //- Read the cell-centre rate data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
