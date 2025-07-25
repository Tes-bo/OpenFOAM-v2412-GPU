/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
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
    Foam::fv::atmPlantCanopyTSource

Group
    grpFvOptionsSources

Description
    Applies sources on temperature, i.e. \c T, to incorporate
    effects of plant canopy for atmospheric boundary layer modelling.

    Corrections applied to:
    \verbatim
      T             | Temperature             [K]
    \endverbatim

    Required fields:
    \verbatim
      T             | Temperature                             [K]
      qPlant        | Tree-height based specific heat flux    [m2/s3]
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmPlantCanopyTSource1
    {
        // Mandatory entries (unmodifiable)
        type                atmPlantCanopyTSource;

        atmPlantCanopyTSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode   all;

            // Optional entries (unmodifiable)
            rho             rho;
            Cp0             1005.0;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                      | Type   | Req'd  | Dflt
      type      | Type name: atmPlantCanopyTSource | word   | yes    | -
      rho       | Name of density field            | word   | no     | rho
      Cp0       | Specific heat capacity [m2/s2/K] | scalar | no     | 1005.0
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    atmPlantCanopyTSource.C
    atmPlantCanopyTSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmPlantCanopyTSource_H
#define fv_atmPlantCanopyTSource_H

#include "cellSetOption.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class atmPlantCanopyTSource Declaration
\*---------------------------------------------------------------------------*/

class atmPlantCanopyTSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Name of density field
        const word rhoName_;

        //- Specific heat capacity
        const dimensionedScalar Cp0_;

        // Fields

            //- Heat flux
            volScalarField qPlant_;


public:

    //- Runtime type information
    TypeName("atmPlantCanopyTSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmPlantCanopyTSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmPlantCanopyTSource(const atmPlantCanopyTSource&) = delete;

        //- No copy assignment
        void operator=(const atmPlantCanopyTSource&) = delete;


    // Member Functions

        //- Add explicit contribution to temperature equation
        //- for incompressible flow computations
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to temperature equation
        //- for compressible flow computations
         virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to temperature equation
        //- for multiphase flow computations
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Read source dictionary (effectively no-op)
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
