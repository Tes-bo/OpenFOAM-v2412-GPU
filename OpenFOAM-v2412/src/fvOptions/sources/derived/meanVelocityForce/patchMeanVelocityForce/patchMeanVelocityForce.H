/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::fv::patchMeanVelocityForce

Description
    Applies the force over the specified patch to maintain
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
    patchMeanVelocityForce1
    {
        // Mandatory entries (unmodifiable)
        type                 patchMeanVelocityForce;
        patch                <patchName>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: patchMeanVelocityForce | word   | yes   | -
      patch     | Name of operand patch             | word   | yes   | -
    \endtable

    The inherited entries are elaborated in:
     - \link meanVelocityForce.H \endlink

See also
  - Foam::fv::meanVelocityForce

SourceFiles
    patchMeanVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef patchMeanVelocityForce_H
#define patchMeanVelocityForce_H

#include "meanVelocityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class patchMeanVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class patchMeanVelocityForce
:
    public meanVelocityForce
{
protected:

    // Protected Data

        //- Name of operand patch
        word patch_;

        //- Patch index
        label patchi_;


    // Protected Member Functions

        //- Calculate and return the magnitude of the mean velocity
        //- averaged over the specified patch
        virtual scalar magUbarAve(const volVectorField& U) const;


public:

    //- Runtime type information
    TypeName("patchMeanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        patchMeanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        patchMeanVelocityForce(const patchMeanVelocityForce&) = delete;

        //- No copy assignment
        void operator=(const patchMeanVelocityForce&) = delete;


    //- Destructor
    ~patchMeanVelocityForce() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
