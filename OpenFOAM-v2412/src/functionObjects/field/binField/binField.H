/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::functionObjects::binField

Group
    grpFieldFunctionObjects

Description
    Calculates binned data, where specified patches are divided into
    segments according to various input bin characteristics, so that
    spatially-localised information can be output for each segment.

    Operands:
    \table
      Operand      | Type                          | Location
      input        | vol\<Type\>Field(s)          <!--
               --> | \<time\>/\<inpField\>s
      output file  | dat                           <!--
               --> | postProcessing/\<FO\>/\<time\>/\<file\>
      output field | -                             | -
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    binField1
    {
        // Mandatory entries
        type                    binField;
        libs                    (fieldFunctionObjects);
        binModel                <word>;
        fields                  (<wordHashSet>);
        patches                 (<wordRes>);
        binData
        {
            // Entries of the chosen binModel
        }

        // Optional entries
        cellZones               (<wordRes>);
        decomposePatchValues    <bool>;

        // Conditional optional entries

            // General coordinate system specification
            coordinateSystem
            {
                type        cartesian;
                origin      (0 0 0);
                rotation
                {
                    type        ...
                    ...
                }
            }

            // Define the centre of rotation
            // with implicit directions e1=(1 0 0) and e3=(0 0 1)
            CofR            (0 0 0);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd | Deflt
      type      | Type name: binField                   | word   | yes  | -
      libs      | Library name: fieldFunctionObjects    | word   | yes  | -
      binModel  | Name of the bin model                 | word   | yes  | -
      fields    | Names of operand fields           | wordHasSet | yes  | -
      patches   | Names of operand patches             | wordRes | yes  | -
      binData   | Entries of the chosen bin model       | dict   | yes  | -
      decomposePatchValues | Flag to output normal and tangential      <!--
                --> components                          | bool   | no   | false
      cellZones | Names of operand cell zones          | wordRes | no   | -
      coordinateSystem | Coordinate system specifier    | dict | cndtnl | -
      CofR    | Centre of rotation                      | vector | cndtnl | -
    \endtable

    Options for the \c binModel entry:
    \verbatim
      singleDirectionUniformBin    | Segments in a single direction
      uniformBin                   | Segments in multiple directions
    \endverbatim

    The inherited entries are elaborated in:
      - \link fvMeshFunctionObject.H \endlink
      - \link writeFile.H \endlink
      - \link coordinateSystem.H \endlink

Note
  - If a \c coordinateSystem entry exists, it is taken in favour of \c CofR.

SourceFiles
    binField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_binField_H
#define Foam_functionObjects_binField_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class binModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class binField Declaration
\*---------------------------------------------------------------------------*/

class binField
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Runtime-selectable bin model
        autoPtr<binModel> binModelPtr_;


public:

    //- Runtime type information
    TypeName("binField");


    // Constructors

        //- Construct from Time and dictionary
        binField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        binField
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        binField(const binField&) = delete;

        //- No copy assignment
        void operator=(const binField&) = delete;


    //- Destructor
    virtual ~binField() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
