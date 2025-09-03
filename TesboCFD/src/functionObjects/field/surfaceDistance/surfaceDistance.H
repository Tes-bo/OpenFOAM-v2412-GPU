/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::functionObjects::surfaceDistance

Group
    grpFieldFunctionObjects

Description
    Computes the distance to the nearest surface from a given geometry.

    Operands:
    \table
      Operand        | Type            | Location
      input          | -               | -
      output file    | -               | -
      output field   | volScalarField  | $FOAM_CASE/\<time\>/surfaceDistance
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    surfaceDistance1
    {
        // Mandatory entries (unmodifiable)
        type            surfaceDistance;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        geometry
        {
            motorBike.obj
            {
                type triSurfaceMesh;
                name motorBike;
            }
        }

        // Optional entries (runtime modifiable)
        calculateCells  true;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: surfaceDistance            | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      geometry     | Surface details                       | dict |  yes  | -
      calculateCells  | Calculate distance from cell       | bool |  no   | true
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::surfaceDistance

SourceFiles
    surfaceDistance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceDistance_H
#define functionObjects_surfaceDistance_H

#include "fvMeshFunctionObject.H"
#include "searchableSurfaces.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class surfaceDistance Declaration
\*---------------------------------------------------------------------------*/

class surfaceDistance
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Switch to calculate distance-to-cells
        Switch doCells_;

        //- Geometry
        autoPtr<searchableSurfaces> geomPtr_;


public:

    //- Runtime type information
    TypeName("surfaceDistance");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceDistance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        surfaceDistance(const surfaceDistance&) = delete;

        //- No copy assignment
        void operator=(const surfaceDistance&) = delete;


    //- Destructor
    virtual ~surfaceDistance() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the interpolated fields
        virtual bool execute();

        //- Write the interpolated fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
