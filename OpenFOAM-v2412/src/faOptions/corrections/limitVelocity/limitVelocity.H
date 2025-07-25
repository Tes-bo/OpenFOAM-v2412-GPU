/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::fa::limitVelocity

Description
    Limits the maximum velocity magnitude to the specified \c max value.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    <faOption>
    {
        // Mandatory entries
        type            limitVelocity;
        active          yes;
        selectionMode   all;
        max             <maxValue>;

        // Optional entries
        U               <UName>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type   | Reqd  | Deflt
      type      | Type name: limitVelocity      | word   | yes   | -
      max       | Maximum velocity limit [m/s]  | scalar | yes   | -
      U         | Name of operand velocity field | word  | no    | U
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link faceSetOption.H \endlink

SourceFiles
    limitVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_limitVelocity_H
#define Foam_fa_limitVelocity_H

#include "faceSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class limitVelocity Declaration
\*---------------------------------------------------------------------------*/

class limitVelocity
:
    public fa::faceSetOption
{
protected:

    // Protected Data

        //- Name of operand velocity field
        word UName_;

        //- Maximum velocity magnitude
        scalar max_;


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
        virtual void correct(areaVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
