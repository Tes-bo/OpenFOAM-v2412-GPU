/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020 OpenCFD Ltd
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
    Foam::fv::atmAmbientTurbSource

Group
    grpFvOptionsSources

Description
    Applies sources on \c k and either \c epsilon or \c omega to prevent them
    droping below a specified ambient value for atmospheric boundary
    layer modelling.  Such adjustment reportedly increases numerical
    stability for very stable atmospheric stability conditions, and prevents
    nonphysical oscillations in regions of low shear at higher altitudes.

    Corrections applied to:
    \verbatim
      k         | Turbulent kinetic energy                  [m2/s2]
    \endverbatim

    Corrections applied to either of the below, if exist:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    Required fields:
    \verbatim
      k             | Turbulent kinetic energy                    [m2/s2]
      epsilon/omega | Dissipation rate OR Spec. dissipation rate  [m2/s3]/[1/s]
    \endverbatim

    References:
    \verbatim
        Background turbulence (tag:RS):
            Rumsey, C. L., & Spalart, P. R. (2009).
            Turbulence model behavior in low Reynolds number
            regions of aerodynamic flowfields.
            AIAA journal, 47(4), 982-993.
            DOI:10.2514/1.39947
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmAmbientTurbSource1
    {
        // Mandatory entries (unmodifiable)
        type                  atmAmbientTurbSource;

        atmAmbientTurbSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Mandatory entries (unmodifiable)
            kAmb              0.0;

            // Optional entries (unmodifiable)
            rho               rho;
            epsilonAmb        0.0;
            omegaAmb          0.0;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                     | Type   | Req'd  | Dflt
      type       | Type name: atmAmbientTurbSource | word   | yes    | -
      kAmb       | Ambient value for \c k          | scalar | yes    | -
      rho        | Name of density field           | word   | no     | rho
      epsilonAmb | Ambient value for \c epsilon    | scalar | no     | 0.0
      omegaAmb   | Ambient value for \c omega      | scalar | no     | 0.0
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    atmAmbientTurbSource.C
    atmAmbientTurbSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmAmbientTurbSource_H
#define fv_atmAmbientTurbSource_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class atmAmbientTurbSource Declaration
\*---------------------------------------------------------------------------*/

class atmAmbientTurbSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Internal flag to determine the working field is epsilon or omega
        Switch isEpsilon_;

        //- Name of density field
        const word rhoName_;

        //- Ambient turbulent kinetic energy value [m2/s2]
        const dimensionedScalar kAmb_;

        //- Ambient turbulent kinetic energy dissipation rate value [m2/s3]
        const dimensionedScalar epsilonAmb_;

        //- Ambient specific dissipation rate value [1/s]
        const dimensionedScalar omegaAmb_;

        //- Required turbulence model coefficients (copied from turb model)
        dimensionedScalar Cmu_;
        dimensionedScalar C2_;


     // Private Member Functions

        //- Apply atmAmbientTurbSource to epsilon field
        template<class AlphaFieldType, class RhoFieldType>
        void atmAmbientTurbSourceEpsilon
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

        //- Apply atmAmbientTurbSource to omega field
        template<class AlphaFieldType, class RhoFieldType>
        void atmAmbientTurbSourceOmega
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

        //- Apply atmAmbientTurbSource to k field
        template<class AlphaFieldType, class RhoFieldType>
        void atmAmbientTurbSourceK
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;


public:

    //- Runtime type information
    TypeName("atmAmbientTurbSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmAmbientTurbSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmAmbientTurbSource(const atmAmbientTurbSource&) = delete;

        //- No copy assignment
        void operator=(const atmAmbientTurbSource&) = delete;


    // Member Functions

        //- Add explicit contribution to epsilon or omega equation
        //- for incompressible flow computations
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for compressible flow computations
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
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

#ifdef NoRepository
    #include "atmAmbientTurbSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
