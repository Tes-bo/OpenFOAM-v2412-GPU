/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::fv::velocityDampingConstraint

Group
    grpFvOptionsConstraints

Description
    Constrain given velocity fields to dampen velocity
    fluctuations exceeding a given value within a specified region.

    This constraint is primarily used to dampen velocity fluctuations in
    the start-up phase of simulations. When the local velocity magnitude
    exceeds the user-supplied maximum value a sink term is activated in
    the affected region to lower the velocity to the limiting value.

    Constraints applied to:
    \verbatim
      U    | Velocity                   [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U    | Velocity                   [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    velocityDampingConstraint1
    {
        // Mandatory entries (unmodifiable)
        type            velocityDampingConstraint;

        // Mandatory entries (runtime modifiable)
        UMax            200;

        // Optional entries (runtime modifiable)
        UNames          (<Uname1> <Uname2> ... <UnameN>);
        C               <scalar>;

        // Conditional optional entries (runtime modifiable)

            // when UNames entry is not present
            U               U;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd  | Dflt
      type      | Type name: velocityDampingConstraint  | word   | yes   | -
      UMax      | Maximum velocity magnitude            | scalar | yes   | -
      UNames    | Names of operand velocity fields      | wordList | no  | -
      U         | Name of operand velocity field        | word   | cndtnl | U
      C         | Model factor                          | scalar | no    | 1
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink

Note
  - When active, this constraint manipulates the system of equations.
    Users should ensure that it is not active when the case is converged
    (steady-state) or during the period of interest (transient) to ensure
    that its presence does not pollute the results.

See also
  - Foam::fv::limitVelocity

SourceFiles
    velocityDampingConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fv_velocityDampingConstraint_H
#define fv_velocityDampingConstraint_H

#include "cellSetOption.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class velocityDampingConstraint Declaration
\*---------------------------------------------------------------------------*/

class velocityDampingConstraint
:
    public fv::cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Protected Data

        //- Maximum velocity magnitude
        scalar UMax_;

        //- Model factor
        scalar C_;


    // Protected Member Functions

        //- Constrain the given velocity fields by a given maximum value
        void addDamping(fvMatrix<vector>& eqn);

        //- Write file header information
        void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("velocityDampingConstraint");


    // Constructors

        //- Construct from components
        velocityDampingConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        velocityDampingConstraint(const velocityDampingConstraint&) = delete;

        //- No copy assignment
        void operator=(const velocityDampingConstraint&) = delete;


    //- Destructor
    virtual ~velocityDampingConstraint() = default;


    // Member Functions

        // Set values directly

            //- Constrain vector matrix
            virtual void constrain(fvMatrix<vector>& eqn, const label fieldi);


        // I-O

            //- Write data
            virtual void writeData(Ostream& os) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
