/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::momentumError

Group
    grpForcesFunctionObjects

Description
    Computes balance terms for the steady momentum equation.

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | -              | -
      output field  | volVectorField | $FOAM_CASE/\<time\>/\<file\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    momentumError1
    {
        // Mandatory entries (unmodifiable)
        type        momentumError;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        p           <pName>;
        U           <UName>;
        phi         <phiName>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd  | Deflt
      type         | Type name: momentumError           | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      p            | Name of pressure field             | word |  no   | p
      U            | Name of velocity field             | word |  no   | U
      phi          | Name of flux field                 | word |  no   | phi
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link zoneSubSet.H \endlink

    Usage by the \c postProcess utility is not available.

Note
  - Optionally the user can specify \c cellZones to create a sub-mesh for the
    \c momentumError calculation.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::momentumError

SourceFiles
    momentumError.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_momentumError_H
#define functionObjects_momentumError_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "zoneSubSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class momentumError Declaration
\*---------------------------------------------------------------------------*/

class momentumError
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        // Read from dictionary

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of flux field
            word phiName_;

            //- Sub-set mesh
            autoPtr<Detail::zoneSubSet> zoneSubSetPtr_;


    // Protected Member Functions

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volVectorField> divDevRhoReff();


public:

    //- Runtime type information
    TypeName("momentumError");


    // Constructors

        //- Construct from Time and dictionary
        momentumError
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        momentumError(const momentumError&) = delete;

        //- No copy assignment
        void operator=(const momentumError&) = delete;


    //- Destructor
    virtual ~momentumError() = default;


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();

        //- Calculate the momentum error
        void calcMomentError();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
