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
    Foam::fv::atmLengthScaleTurbSource

Group
    grpFvOptionsSources

Description
    Applies sources on either \c epsilon or \c omega to correct
    mixing-length scale estimations for atmospheric boundary layer modelling.

    Corrections applied to either of the below, if exist:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    Required fields, either of the below:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    References:
    \verbatim
        Mixing-length scale limiter for epsilon (tag:AC):
            Apsley, D. D., & Castro, I. P. (1997).
            A limited-length-scale k-ε model for the neutral and
            stably-stratified atmospheric boundary layer.
            Boundary-layer meteorology, 83(1), 75-98.
            DOI:10.1023/A:1000252210512

        Mixing-length scale limiter for omega (tag:L):
            Langner, J. (2016).
            Implementierung und validierung von RANS-modellen der
            thermisch geschichteten, atmosphärischen grenzschicht.
            Masterarbeit zum thema, Technische Universität Berlin.

        Mixing-length scale estimation (tag:P):
            Pope, S. B. (2000).
            Turbulent flows.
            Cambridge, UK: Cambridge Univ. Press
            DOI:10.1017/CBO9780511840531
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmLengthScaleTurbSource1
    {
        // Mandatory entries (unmodifiable)
        type                 atmLengthScaleTurbSource;

        atmLengthScaleTurbSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Optional entries (unmodifiable)
            rho              rho;
            Lmax             41.575;
            n                3.0;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type   | Req'd  | Dflt
      type      | Type name: atmLengthScaleTurbSource | word   | yes    | -
      rho       | Name of density field               | word   | no     | rho
      Lmax      | Maximum mixing-length scale [m]     | scalar | no     | 41.575
      n         | Mixing-length scale exponent        | scalar | no     | 3.0
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

See also
    - cellSetOption.H

SourceFiles
    atmLengthScaleTurbSource.C
    atmLengthScaleTurbSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmLengthScaleTurbSource_H
#define fv_atmLengthScaleTurbSource_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class atmLengthScaleTurbSource Declaration
\*---------------------------------------------------------------------------*/

class atmLengthScaleTurbSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Internal flag to determine the working field is epsilon or omega
        bool isEpsilon_;

        //- Name of density field
        const word rhoName_;

        //- Maximum mixing-length scale [m]
        const dimensionedScalar Lmax_;

        //- Mixing-length scale exponent
        const dimensionedScalar n_;

        //- Required turbulence model coefficients (copied from turb model)
        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar C3_;


     // Private Member Functions

        //- Return mixing-length scale modifier for epsilon-based models
        tmp<volScalarField::Internal> calcC1Star
        (
            const volScalarField::Internal& k,
            const volScalarField::Internal& epsilon
        ) const;

        //- Return mixing-length scale modifier for omega-based models
        tmp<volScalarField::Internal> calcGammaStar
        (
            const volScalarField::Internal& k,
            const volScalarField::Internal& omega,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;

        //- Apply atmLengthScaleTurbSource to epsilon
        template<class AlphaFieldType, class RhoFieldType>
        void atmLengthScaleTurbSourceEpsilon
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

        //- Apply atmLengthScaleTurbSource to omega
        template<class AlphaFieldType, class RhoFieldType>
        void atmLengthScaleTurbSourceOmega
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;


public:

    //- Runtime type information
    TypeName("atmLengthScaleTurbSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmLengthScaleTurbSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmLengthScaleTurbSource(const atmLengthScaleTurbSource&) = delete;

        //- No copy assignment
        void operator=(const atmLengthScaleTurbSource&) = delete;


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
    #include "atmLengthScaleTurbSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
