/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::fv::limitTurbulenceViscosity

Group
    grpFvOptionsCorrections

Description
    Corrects turbulence viscosity field (e.g. \c nut) within a specified
    region by applying a maximum limit, set according to a coefficient
    multiplied by the laminar viscosity:

    \f[
        \nu_{t,max} = c \nu
    \f]

    Corrections applied to:
    \verbatim
      nut  | Turbulence vicosity                   [m2/s2]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    limitTurbulenceViscosity1
    {
        // Mandatory entries (unmodifiable)
        type            limitTurbulenceViscosity;

        // Optional entries (runtime modifiable)
        nut             nut;
        c               1e5;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                | Type   | Reqd  | Dflt
      type      | Type name: limitTurbulenceViscosity | word  | yes   | -
      nut       | Name of nut field          | word   | no    | nut
      c         | Limiting coefficient       | scalar | yes   | 1e5
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink

SourceFiles
    limitTurbulenceViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef fv_limitTurbulenceViscosity_H
#define fv_limitTurbulenceViscosity_H

#include "cellSetOption.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class limitTurbulenceViscosity Declaration
\*---------------------------------------------------------------------------*/

class limitTurbulenceViscosity
:
    public fv::cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Protected Data

        //- Name of turbulence viscosity field
        word nutName_;

        //- Limiting coefficient []
        scalar c_;


    // Protected Member Functions

        //- Write file header information
        void writeFileHeader(Ostream& os);

        //- Return laminar viscosity
        tmp<volScalarField> nu() const;


public:

    //- Runtime type information
    TypeName("limitTurbulenceViscosity");


    // Constructors

        //- Construct from components
        limitTurbulenceViscosity
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        limitTurbulenceViscosity(const limitTurbulenceViscosity&) = delete;

        //- No copy assignment
        void operator=(const limitTurbulenceViscosity&) = delete;


    //- Destructor
    virtual ~limitTurbulenceViscosity() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the energy field
        virtual void correct(volScalarField& nut);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
