/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::fv::explicitPorositySource

Group
    grpFvOptionsSources

Description
    Applies an explicit porosity source to
    the momentum equation within a specified region.

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
    explicitPorositySource1
    {
        // Mandatory entries (unmodifiable)
        type                 <porosityModelName>;

        <porosityModelName>Coeffs
        {
            // Optional entries (runtime modifiable)
            UNames          (<Uname1> <Uname2> ... <UnameN>);

            // Conditional optional entries (runtime modifiable)

                // when the entry "UNames" is not present
                U           U;

            // Mandatory/Optional (inherited) entries
            ...
        }

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: \<porosityModelName\>  | word   | yes   | -
      UNames    | Names of operand velocity fields  | wordList | no  | -
      U         | Name of operand velocity field    | word   | no    | U
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

    Example usage, here employing the Darcy-Forchheimer model:
    \verbatim
    explicitPorositySourceCoeffs
    {
        type            DarcyForchheimer;

        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);

            coordinateSystem
            {
                origin  (0 0 0);
                e1      (0.70710678 0.70710678 0);
                e2      (0 0 1);
            }
        }
    }
    \endverbatim

Note
  - The porous region must be selected as a \c cellZone.

See also
  - Foam::porosityModel

SourceFiles
    explicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_explicitPorositySource_H
#define Foam_explicitPorositySource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class explicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class explicitPorositySource
:
    public fv::cellSetOption
{
protected:

    // Protected Data

        //- Run-time selectable porosity model
        autoPtr<porosityModel> porosityPtr_;


public:

    //- Runtime type information
    TypeName("explicitPorositySource");


    // Constructors

        //- Construct from components
        explicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        explicitPorositySource(const explicitPorositySource&) = delete;

        //- No copy assignment
        void operator=(const explicitPorositySource&) = delete;


    //- Destructor
    virtual ~explicitPorositySource() = default;


    // Member Functions

        //- Access to the porosityModel
        const porosityModel& model() const noexcept
        {
            return *porosityPtr_;
        }


        //- Add implicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add implicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add implicit contribution to phase momentum equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
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
