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
    Foam::functionObjects::wallShearStress

Group
    grpForcesFunctionObjects

Description
    Computes the wall-shear stress at selected wall patches.

    \f[
        \vec \tau = \vec R \cdot \vec n
    \f]

    where
    \vartable
        \vec \tau   | Wall shear stress vector [m^2/s^2]
        \vec R | Shear-stress symmetric tensor (retrieved from turbulence model)
        \vec n | Patch normal vector (into the domain)
    \endvartable

    Operands:
    \table
      Operand        | Type              | Location
      input          | -                 | -
      output file    | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<field\>
      output field   | volVectorField (only boundaryField) <!--
                 --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallShearStress1
    {
        // Mandatory entries (unmodifiable)
        type        wallShearStress;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");

        // Optional entries
        writeFields     true;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: wallShearStress            | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patches    | Names of operand patches   | wordRes | no | all wall patches
      writeFields | Flag to write shear stress field       | bool |  no   | true
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func wallShearStress
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
    - ExtendedCodeGuide::functionObjects::field::wallShearStress

SourceFiles
    wallShearStress.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_wallShearStress_H
#define Foam_functionObjects_wallShearStress_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallShearStress Declaration
\*---------------------------------------------------------------------------*/

class wallShearStress
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Write the shear stress field ?
        bool writeFields_;


protected:

    // Protected Data

        //- Wall patches to process (optionally filtered by name)
        labelList patchIDs_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the shear-stress
        void calcShearStress
        (
            const volSymmTensorField& Reff,
            volVectorField& shearStress
        );


public:

    //- Runtime type information
    TypeName("wallShearStress");


    // Constructors

        //- Construct from Time and dictionary
        wallShearStress
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        wallShearStress(const wallShearStress&) = delete;

        //- No copy assignment
        void operator=(const wallShearStress&) = delete;


    //- Destructor
    virtual ~wallShearStress() = default;


    // Member Functions

        //- Read the wallShearStress data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Report min/max and log to file,
        //- write the wall shear-stress volume field.
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
