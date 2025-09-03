/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::functionObjects::momentum

Group
    grpFieldFunctionObjects

Description
    Computes linear/angular momentum, reporting integral values
    and optionally writing the fields.

    Operands:
    \table
      Operand       | Type | Location
      input         | -    | -
      output file   | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field  | -    | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    momentum1
    {
        // Mandatory entries (unmodifiable)
        type            momentum;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        regionType      all;
        writeMomentum   yes;
        writePosition   yes;
        writeVelocity   yes;
        p               p;
        U               U;
        rho             rho;
        rhoRef          1.0;

        cylindrical     true;
        origin          (0 0 0);
        e1              (1 0 0);
        e3              (0 0 1);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: momentum                | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      regionType   | Selection type: all/cellSet/cellZone | word   | no  | all
      writeMomentum | Write (linear, angular) momentum  fields | bool | no | no
      writePosition | Write angular position component fields | bool | no  | no
      writeVelocity | Write angular velocity fields           | bool | no  | no
      p            | Pressure field name                | word | no  | p
      U            | Velocity field name                | word | no  | U
      rho          | Density field name                 | word | no  | rho
      rhoRef       | Reference density (incompressible) | scalar | no  | 1.0
      cylindrical  | Use cylindrical coordinates        | bool | no  | no
      origin  | Origin for cylindrical coordinates   | vector | conditional  | -
      name    | Name of cellSet/cellZone if required | word   | conditional  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

Note
  - For incompressible cases, the value of \c rhoRef is used.
  - When specifying the cylindrical coordinate system, the rotation
    can be specified directly with e1/e2/e3 axes, or via a \c rotation
    sub-dictionary
    For example,
    \verbatim
        origin      (0 0 0);
        rotation
        {
            type    cylindrical;
            axis    (0 0 1);
        }
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::volRegion
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::momentum

SourceFiles
    momentum.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_momentum_H
#define functionObjects_momentum_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "cylindricalCS.H"
#include "volFieldsFwd.H"
#include "volRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dimensionSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class momentum Declaration
\*---------------------------------------------------------------------------*/

class momentum
:
    public fvMeshFunctionObject,
    public volRegion,
    public writeFile
{
    // Private Member Functions

        //- Remove calculated fields from the registry
        void purgeFields();

        //- Calculate the fields and integral values
        void calc();

        //- Allocate a new zero geometric field
        template<class GeoField>
        autoPtr<GeoField> newField
        (
            const word& baseName,
            const dimensionSet& dims,
            bool registerObject=true
        ) const;


protected:

    // Protected Data

        //- Integral (linear) momentum
        vector sumMomentum_;

        //- Integral angular momentum
        vector sumAngularMom_;


    // Read from dictionary

        //- The velocity field name (optional)
        word UName_;

        //- The pressure field name (optional)
        //  Only used to determine incompressible/compressible
        word pName_;

        //- The density field name (optional)
        word rhoName_;

        //- Reference density (for incompressible)
        scalar rhoRef_;

        //- Coordinate system for evaluating angular momentum
        coordSystem::cylindrical csys_;

        //- Are we using the cylindrical coordinate system?
        bool hasCsys_;

        //- Write fields flag
        bool writeMomentum_;

        //- Write fields flag
        bool writeVelocity_;

        //- Write fields flag
        bool writePosition_;

        //- Initialised flag
        bool initialised_;


    // Protected Member Functions

        //- Initialise the fields
        void initialise();

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Write momentum data
        void writeValues(Ostream& os);


public:

    //- Runtime type information
    TypeName("momentum");


    // Constructors

        //- Construct from Time and dictionary
        momentum
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        momentum
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        momentum(const momentum&) = delete;

        //- No copy assignment
        void operator=(const momentum&) = delete;


    //- Destructor
    virtual ~momentum() = default;


    // Member Functions

        //- Read the momentum data
        virtual bool read(const dictionary&);

        //- Calculate and report the integral momentum
        virtual bool execute();

        //- Write the momentum, possibly angular momentum and velocity
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
