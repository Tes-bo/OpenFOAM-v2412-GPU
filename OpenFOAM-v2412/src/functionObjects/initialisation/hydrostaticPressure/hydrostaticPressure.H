/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::functionObjects::hydrostaticPressure

Group
    grpInitialisationFunctionObjects

Description
    Calculates and outputs the pressure fields p_rgh and ph_rgh

Usage
    Example of function object specification to calculate hydrostatic pressure:
    \verbatim
    hydrostaticPressure1
    {
        type        hydrostaticPressure;
        libs        (initialisationFunctionObjects);
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        log          | Log to standard output  | no          | yes
        p_rgh        | Name of p_rgh field     | no          | p_rgh
        ph_rgh       | Name of ph_rgh field    | no          | ph_rgh
        pRef         | Name of pressure ref field | no       | pRef
        pRefValue    | Value of pressure ref | when pRef set to pInf | 0
        rho          | Name of density field   | no          | rho
        U            | Name of velocity field  | no          | U
        gh           | Name of gravity*height volume field  | no | gh
        ghf          | Name of gravity*height surface field | no | ghf
        nCorrectors  | Number of correctors when solving ph_rgh | no | 5
    \endtable

Note
    Calculates the hydrostatic pressure on construction/re-initialisation;
    the execute and write functions are not used.

SourceFiles
    hydrostaticPressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_hydrostaticPressure_H
#define functionObjects_hydrostaticPressure_H

#include "fvMeshFunctionObject.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class hydrostaticPressure Declaration
\*---------------------------------------------------------------------------*/

class hydrostaticPressure
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Name of p - rho*g*h field, default is "p_rgh"
        word p_rghName_;

        //- Name of p_hydrostatic - rho*g*h field, default is "ph_rgh"
        word ph_rghName_;

        //- Name of uniform pressure reference field, default is "pRef"
        //  Set to "none" to set to zero or pInf to set according to pRefValue
        word pRefName_;

        //- Reference pressure if pRefName is set to "pInf"
        scalar pRefValue_;

        //- Name of density field, default is "rho"
        word rhoName_;

        //- Name of velocity field, default is "ph_rgh"
        word UName_;

        //- Name of g*h volume field, default is "gh"
        word ghName_;

        //- Name of g*h surface field, default is "ghf"
        word ghfName_;

        //- Number of correctors when solving for ph_rgh
        label nCorrectors_;


    // Protected Member Functions

        //- Helper function to return the reference pressure
        //  Returns zero if pRefName = none; otherwise returns the value
        //  from the mesh database
        dimensionedScalar pRef() const;

        //- Calculate the fields and write
        void calculateAndWrite();


public:

    //- Runtime type information
    TypeName("hydrostaticPressure");


    // Constructors

        //- Construct from Time and dictionary
        hydrostaticPressure
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~hydrostaticPressure() = default;


    // Member Functions

        //- Read the hydrostaticPressure data
        virtual bool read(const dictionary& dict);

        //- Calculate the p_rgh field
        virtual bool execute();

        //- Write the p_rgh and derived fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
