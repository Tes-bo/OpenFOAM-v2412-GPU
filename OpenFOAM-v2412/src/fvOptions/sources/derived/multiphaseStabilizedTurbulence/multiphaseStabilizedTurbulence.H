/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd
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
    Foam::fv::multiphaseStabilizedTurbulence

Group
    grpFvOptionsSources

Description
    Applies corrections to the turbulent kinetic energy equation (i.e. \c k)
    and turbulent viscosity field (i.e. \c nut) for incompressible multiphase
    flow cases.

    Turbulent kinetic energy is over-predicted in VOF solvers at the phase
    interface and throughout the water column in nearly-potential flow regions
    beneath surface waves.

    References:
    \verbatim
        Buoyancy source term in turbulent kinetic energy equation:
            Devolder, B., Rauwoens, P., and Troch, P. (2017).
            Application of a buoyancy-modified k-w SST turbulence model to
            simulate wave run-up around a monopile subjected to regular waves
            using OpenFOAM.
            Coastal Engineering, 125, 81-94.
            DOI:10.1016/j.coastaleng.2017.04.004

        Correction to turbulent viscosity:
            Larsen, B.E. and Fuhrman, D.R. (2018).
            On the over-production of turbulence beneath surface waves in
            Reynolds-averaged Navier-Stokes models
            J. Fluid Mech, 853, 419-460.
            DOI:10.1017/jfm.2018.577
    \endverbatim

    Sources applied to:
    \verbatim
      k         | Turbulent kinetic energy                   [m2/s2]
      nut       | Turbulent viscosity                        [m/s2]
    \endverbatim

    Required fields:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate  [m2/s3]
      k         | Turbulent kinetic energy                   [m2/s2]
      nut       | Turbulent viscosity                        [m/s2]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    multiphaseStabilizedTurbulence1
    {
        // Mandatory entries (unmodifiable)
        type                 multiphaseStabilizedTurbulence;

        // Optional entries (unmodifiable)
        rho             <rhoName>;
        Cmu             0.09;  // from k-epsilon model
        C               1.51;  // from k-omega model
        lambda2         0.1;   // A value of 0 sets the nut correction to 0
        alpha           1.36;  // 1/Prt

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd | Dflt
      type      | Type name: multiphaseStabilizedTurbulence  | word | yes | -
      rho       | Name of operand density field     | word   | yes  | rho
      Cmu       | Model coefficient                 | scalar | no   | 0.09
      C         | Model coefficient                 | scalar | no   | 1.51
      lambda2   | Ratio of nut correction           | scalar | no   | 0.1
      alpha     | Inverse of turbulent Prandtl number | scalar | no | 1.36
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

Note
  - The model C coefficient for the k-epsilon model equates to C2/C1 = 1.33;
    the (default) value of 1.51 comes from the k-omega model and is more
    conservative.
  - Only applicable to epsilon-based turbulence closure models.

SourceFiles
    multiphaseStabilizedTurbulence.C

\*---------------------------------------------------------------------------*/

#ifndef fv_multiphaseStabilizedTurbulence_H
#define fv_multiphaseStabilizedTurbulence_H

#include "fvOption.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class multiphaseStabilizedTurbulence Declaration
\*---------------------------------------------------------------------------*/

class multiphaseStabilizedTurbulence
:
    public fv::option
{
    // Private Data

        //- Name of density field
        const word rhoName_;


        // Model coefficients

            //- k-epsilon model Cmu coefficient
            dimensionedScalar Cmu_;

            //- Model coefficient
            //  For k-epsilon models this equates to C2/C1 = 1.33 and for
            //  k-omega = 1.51.  Here the default is the more conservative 1.51
            dimensionedScalar C_;

            //- lambda2 coefficient; default = 0.1
            dimensionedScalar lambda2_;

            //- Buoyancy coefficient
            dimensionedScalar alpha_;


public:

    //- Runtime type information
    TypeName("multiphaseStabilizedTurbulence");


    // Constructors

        //- Construct from explicit source name and mesh
        multiphaseStabilizedTurbulence
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        multiphaseStabilizedTurbulence
        (
            const multiphaseStabilizedTurbulence&
        ) = delete;

        //- No copy assignment
        void operator=(const multiphaseStabilizedTurbulence&) = delete;


    //- Destructor
    virtual ~multiphaseStabilizedTurbulence() = default;


    // Member Functions

        //- Add explicit contribution to compressible k equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to incompressible k equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Correct the turbulent viscosity
        virtual void correct(volScalarField& field);

        //- Read source dictionary
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
