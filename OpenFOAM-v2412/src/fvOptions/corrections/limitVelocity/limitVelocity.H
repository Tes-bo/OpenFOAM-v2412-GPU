/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::fv::limitVelocity

Description
    Corrects velocity field (i.e. \c T) within a specified
    region by applying a given maximum velocity magnitude limit.

    Corrections applied to:
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
    limitVelocity1
    {
        // Mandatory entries (unmodifiable)
        type            limitVelocity;

        // Mandatory entries (runtime modifiable)
        max             <maxValue>;

        // Optional entries (runtime modifiable)
        U               <Uname>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type   | Reqd  | Dflt
      type      | Type name: limitVelocity      | word   | yes   | -
      max       | Maximum velocity limit [m/s]  | scalar | yes   | -
      U         | Name of operand velocity field | word  | no    | U
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink

See also
  - Foam::fv::velocityDampingConstraint

SourceFiles
    limitVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef fv_limitVelocity_H
#define fv_limitVelocity_H

#include "cellSetOption.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class limitVelocity Declaration
\*---------------------------------------------------------------------------*/

class limitVelocity
:
    public fv::cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Protected Data

        //- Name of operand velocity field
        word UName_;

        //- Maximum velocity magnitude
        scalar max_;


    // Protected Member Functions

        //- Write file header information
        void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("limitVelocity");


    // Constructors

        //- Construct from components
        limitVelocity
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        limitVelocity(const limitVelocity&) = delete;

        //- No copy assignment
        void operator=(const limitVelocity&) = delete;


    //- Destructor
    virtual ~limitVelocity() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the velocity field
        virtual void correct(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
