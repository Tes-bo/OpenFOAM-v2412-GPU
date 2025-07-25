/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::fv::acousticDampingSource

Group
    grpFvOptionsSources

Description
    Applies sources on velocity (i.e. \c U)
    within a specified region to enable acoustic damping.

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
    acousticDampingSource1
    {
        // Mandatory entries (unmodifiable)
        type                 acousticDampingSource;

        // Mandatory entries (runtime modifiable)
        frequency        3000;
        URef             UMean;
        origin           (-1.25 0 0);
        radius1          1.2;
        radius2          1.65;

        // Optional entries (runtime modifiable)
        w                20;

        // Conditional optional entries (runtime modifiable)
        // Either of the below

            // Option-1
            UNames           (<UName1> <UName2> ... <UNameN>);

            // Option-2
            U                U;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: acousticDampingSource  | word   | yes   | -
      frequency | Frequency value [Hz]              | scalar | yes   | -
      URef      | Name of reference velocity field  | word   | yes   | -
      origin    | Sphere centre location or damping | vector | yes   | -
      radius1   | Inner radius at which to start damping | scalar | yes | -
      radius2 | Outer radius beyond which damping is applied | scalar | yes | -
      w         | Stencil width                     | label  | no    | 20
      UNames    | Names of operand velocity fields  | wordList | cndtnl | -
      U         | Name of operand velocity field    | word   | cndtnl   | U
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

SourceFiles
    acousticDampingSource.C

\*---------------------------------------------------------------------------*/

#ifndef acousticDampingSource_H
#define acousticDampingSource_H

#include "cellSetOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class acousticDampingSource Declaration
\*---------------------------------------------------------------------------*/

class acousticDampingSource
:
    public fv::cellSetOption
{
protected:

    // Protected Data

        //- Blending factor [-]
        volScalarField blendFactor_;

        //- Frequency [Hz]
        dimensionedScalar frequency_;

        // Sphere centre location or damping
        point x0_;

        // Inner radius at which to start damping
        scalar r1_;

        // Outer radius beyond which damping is applied
        scalar r2_;

        //- Name of reference velocity field
        word URefName_;

        //- Stencil width
        label w_;


    // Protected Member Functions

        //- Helper function to set the blending factor
        void setBlendingFactor();


public:

    //- Runtime type information
    TypeName("acousticDampingSource");


    // Constructors

        //- Construct from components
        acousticDampingSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        acousticDampingSource(const acousticDampingSource&) = delete;

        //- No copy assignment
        void operator=(const acousticDampingSource&) = delete;


    //- Destructor
    virtual ~acousticDampingSource() = default;


    // Member Functions

        //- Add implicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        //- Add implicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        );

        //- Add implicit contribution to phase momentum equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldI
        );


        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
