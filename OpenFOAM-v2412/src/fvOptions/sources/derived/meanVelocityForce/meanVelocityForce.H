/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::meanVelocityForce

Group
    grpFvOptionsSources

Description
    Applies the force within a specified region to maintain
    the specified mean velocity for incompressible flows.

    Sources applied to:
    \verbatim
      U         | Velocity                                   [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity                                   [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    meanVelocityForce1
    {
        // Mandatory entries (unmodifiable)
        type                 meanVelocityForce;

        // Mandatory entries (unmodifiable)
        Ubar             (1 0 0);
        fields           (<fieldName>);

        // Optional entries (unmodifiable)
        relaxation       1.0;

        // Conditional mandatory entries (unmodifiable)

            // when <timePath>/uniform/<name>Properties file exists
            gradient    <vectorField>; // reading from the aforementioned file

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: meanVelocityForce      | word   | yes   | -
      Ubar      | Desired mean velocity             | vector | yes   | -
      fields    | Name of operand velocity field    | word   | yes   | -
      relaxation | Relaxation factor                | scalar | no    | 1
      gradient  | Initial pressure gradient field   | vectorField | cndtnl | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

Note
  - Currently only handles kinematic pressure (incompressible solvers).

SourceFiles
    meanVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef meanVelocityForce_H
#define meanVelocityForce_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class meanVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class meanVelocityForce
:
    public fv::cellSetOption
{
protected:

    // Protected Data

        //- Desired mean velocity
        vector Ubar_;

        //- Pressure gradient before correction
        scalar gradP0_;

        //- Change in pressure gradient
        scalar dGradP_;

        //- Flow direction
        vector flowDir_;

        //- Relaxation factor
        scalar relaxation_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> rAPtr_;


    // Protected Member Functions

        //- Calculate and return the magnitude of the mean velocity
        //- averaged over the selected cellSet
        virtual scalar magUbarAve(const volVectorField& U) const;

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("meanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        meanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        meanVelocityForce(const meanVelocityForce&) = delete;

        //- No copy assignment
        void operator=(const meanVelocityForce&) = delete;


    //- Destructor
    virtual ~meanVelocityForce() = default;


    // Member Functions

        // Evaluation

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
