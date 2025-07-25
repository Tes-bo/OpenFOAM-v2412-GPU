/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::yPlus

Group
    grpFieldFunctionObjects

Description
    Computes the near wall \f$y_1^+\f$ for turbulence models.

    Operands:
    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    yPlus1
    {
        // Mandatory entries (unmodifiable)
        type        yPlus;
        libs        (fieldFunctionObjects);

        // Optional entries
        useWallFunction true;
        writeFields     true;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                        | Type | Req'd | Dflt
      type       | Type name: yPlus                   | word |  yes  | -
      libs       | Library name: fieldFunctionObjects | word |  yes  | -
      useWallFunction | Flag to use the local expressions of <!--
                      --> the selected nut wall function to compute yPlus, <!--
                      --> otherwise directly compute yPlus from flow field <!--
                      -->                             | bool | no    | true
      writeFields | Flag to write yPlus field         | bool |  no   | true
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func yPlus
    \endverbatim

Note
    The \c writeFields flag currently defaults to \c true
    (for consistency with previous versions that did not have that flag)
    but this is subject to change in the near future for consistency
    with other function objects.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::yPlus

SourceFiles
    yPlus.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_yPlus_H
#define Foam_functionObjects_yPlus_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class yPlus Declaration
\*---------------------------------------------------------------------------*/

class yPlus
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Use the wall function expressions to compute yPlus ?
        bool useWallFunction_;

        //- Write the yPlus field ?
        bool writeFields_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("yPlus");


    // Constructors

        //- Construct from Time and dictionary
        yPlus
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        yPlus(const yPlus&) = delete;

        //- No copy assignment
        void operator=(const yPlus&) = delete;


    //- Destructor
    virtual ~yPlus() = default;


    // Member Functions

        //- Read the yPlus data
        virtual bool read(const dictionary&);

        //- Calculate the yPlus field
        virtual bool execute();

        //- Report min/max/avg (per patch) and log to file,
        //- write the yPlus volume field.
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
