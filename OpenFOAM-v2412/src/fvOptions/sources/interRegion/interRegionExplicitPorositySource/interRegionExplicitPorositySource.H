/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::fv::interRegionExplicitPorositySource

Group
    grpFvOptionsSources

Description
    Applies inter-region explicit porosity source.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    interRegionExplicitPorositySource1
    {
        // Mandatory entries (unmodifiable)
        type                 interRegionExplicitPorositySource;

        interRegionExplicitPorositySourceCoeffs
        {
            // Mandatory entries (unmodifiable)
            type     <porosityModel>;

            // Optional entries (runtime modifiable)
            U        U;
            mu       thermo:mu;
        }

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                     | Type   | Reqd  | Dflt
      type  | Type name: interRegionExplicitPorositySource | word  | yes | -
      type  | Name of porosity model              | word   | yes   | -
      U     | Name of operand velocity field      | word   | no    | U
      mu    | Name of operand dynamic viscosity field | word | no  | thermo:mu
    \endtable

    The inherited entries are elaborated in:
     - \link interRegionOption.H \endlink

    Sources described by, for example using the DarcyForchheimer model:

    \verbatim
    interRegionExplicitPorositySourceCoeffs
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
  - The porous region must be selected as a cellZone.

SourceFiles
    interRegionExplicitPorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionExplicitPorositySource_H
#define interRegionExplicitPorositySource_H

#include "interRegionOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class porosityModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class interRegionExplicitPorositySource Declaration
\*---------------------------------------------------------------------------*/

class interRegionExplicitPorositySource
:
    public interRegionOption
{
protected:

    // Protected Data

        //- Run-time selectable porosity model
        autoPtr<porosityModel> porosityPtr_;

        //- First iteration
        bool firstIter_;

        //- Name of operand velocity field
        word UName_;

        //- Name of operand dynamic viscosity field (compressible case only)
        word muName_;


    // Protected Member Functions

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("interRegionExplicitPorositySource");


    // Constructors

        //- Construct from components
        interRegionExplicitPorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        interRegionExplicitPorositySource
        (
            const interRegionExplicitPorositySource&
        ) = delete;

        //- No copy assignment
        void operator=(const interRegionExplicitPorositySource&) = delete;


    //- Destructor
    virtual ~interRegionExplicitPorositySource() = default;


    // Member Functions

        // Add explicit and implicit contributions

            //- Vector
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // Add explicit and implicit contributions to compressible equation

            //- Vector
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
