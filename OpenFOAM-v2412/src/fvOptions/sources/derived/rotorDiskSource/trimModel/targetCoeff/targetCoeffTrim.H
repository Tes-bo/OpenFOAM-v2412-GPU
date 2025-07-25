/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::targetCoeffTrim

Description
    Trim model where the operating characteristics of rotor
    (e.g. blade pitch angle) can vary to reach a specified
    target thrust and torque.

    Solves:

    \f[
        c^{old} + J.d(\theta) = c^{target}
    \f]

    where
    \vartable
      n      | Time level
      c      | Coefficient vector (thrust force, roll moment, pitch moment)
      \theta | Pitch angle vector (collective, roll, pitch)
      J      | Jacobian [3x3] matrix
    \endvartable

    The trimmed pitch angles are found via solving the above with a
    Newton-Raphson iterative method.  The solver tolerance can be user-input,
    using the 'tol' entry.

    If coefficients are requested (useCoeffs = true),
    the force and moments are normalised using:

    \f[
        c = \frac{F}{\alpha \pi \rho \omega^2 R^4}
    \f]

    and

    \f[
        c = \frac{M}{\alpha \pi \rho \omega^2 R^5}
    \f]

    where
    \vartable
      F      | Force
      M      | Moment
      \alpha | User-input conversion coefficient (default = 1)
      \rho   | Fluid desity
      \omega | Rotor angular velocity
      \pi    | Mathematical pi
      R      | Rotor radius
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSource1
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (runtime modifiable)
        trimModel    targetForce;

        targetForceCoeffs
        {
            // Conditional mandatory entries (runtime modifiable)

                // when trimModel=targetForce
                target
                {
                    // Mandatory entries (runtime modifiable)
                    thrust       0.5;
                    pitch        0.5;
                    roll         0.5;

                    // Optional entries (runtime modifiable)
                    useCoeffs    true;
                }

                pitchAngles
                {
                    theta0Ini    0.1;
                    theta1cIni   0.1;
                    theta1sIni   0.1;
                }

                calcFrequency    1;

                // Optional entries (runtime modifiable)
                nIter            50;
                tol              1e-8;
                relax            1;
                dTheta           0.1;
                alpha            1.0;
        }
    }

    where the entries mean:
    \table
      Property       | Description                   | Type  | Reqd | Dflt
      calcFrequency  | Number of iterations between calls to 'correct' <!--
                     -->                             | label | yes  | -
      useCoeffs      | Flag to indicate whether to solve coeffs <!--
                     --> (true) or forces (false)    | bool  | no   | true
      thrust         | Target thrust (coefficient)   | scalar | yes | -
      pitch          | Target pitch (coefficient)    | scalar | yes | -
      roll           | Target roll moment (coefficient)   | scalar | yes | -
      theta0Ini  | Initial pitch angle [deg]         | scalar | yes | -
      theta1cIni | Initial lateral pitch angle (cos coeff) [deg] <!--
                 -->                                 | scalar | yes | -
      theta1sIni | Initial longitudinal pitch angle (sin coeff) [deg] <!--
                 -->                                 | scalar | yes | -
      nIter      | Maximum number of iterations in trim routine <!--
                 -->                                 | label  | no  | 50
      tol        | Convergence tolerance in trim routine | scalar | no | 1e-8
      relax      | Relaxation factor in trim routine | scalar | no  | 1
      dTheta     | Perturbation angle used to determine Jacobian [deg] <!--
                 -->                                 | scalar | no  | 0.1
      alpha          | Coefficient to allow for conversion between US <!--
                     --> and EU definitions          | scalar | no  | 1
    \endtable

See also
  - Foam::fv::rotorDiskSource
  - Foam::trimModel
  - Foam::fixedTrim

SourceFiles
    targetCoeffTrim.C

\*---------------------------------------------------------------------------*/

#ifndef targetCoeffTrim_H
#define targetCoeffTrim_H

#include "trimModel.H"
#include "tensor.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class targetCoeffTrim Declaration
\*---------------------------------------------------------------------------*/

class targetCoeffTrim
:
    public trimModel
{
protected:

    // Protected Data

        //- Number of iterations between calls to 'correct'
        label calcFrequency_;

        //- Flag to indicate whether to solve coeffs (true) or forces (false)
        bool useCoeffs_;

        //- Target coefficient vector (thrust force, roll moment, pitch moment)
        vector target_;

        //- Pitch angles (collective, roll, pitch) [rad]
        vector theta_;

        //- Maximum number of iterations in trim routine
        label nIter_;

        //- Convergence tolerance
        scalar tol_;

        //- Under-relaxation coefficient
        scalar relax_;

        //- Perturbation angle used to determine jacobian
        scalar dTheta_;

        //- Coefficient to allow for conversion between US and EU definitions
        scalar alpha_;


    // Protected Member Functions

        //- Calculate the rotor force and moment coefficients vector
        template<class RhoFieldType>
        vector calcCoeffs
        (
            const RhoFieldType& rho,
            const vectorField& U,
            const scalarField& alphag,
            vectorField& force
        ) const;

        //- Correct the model
        template<class RhoFieldType>
        void correctTrim
        (
            const RhoFieldType& rho,
            const vectorField& U,
            vectorField& force
        );


public:

    //- Run-time type information
    TypeName("targetCoeffTrim");

    // Constructors

        //- Constructor from rotor and dictionary
        targetCoeffTrim
        (
            const fv::rotorDiskSource& rotor,
            const dictionary& dict
        );

        //- No copy construct
        targetCoeffTrim(const targetCoeffTrim&) = delete;

        //- No copy assignment
        void operator=(const targetCoeffTrim&) = delete;


    //- Destructor
    virtual ~targetCoeffTrim() = default;


    // Member Functions

        //- Read
        void read(const dictionary& dict);

        //- Return the geometric angle of attack [rad]
        virtual tmp<scalarField> thetag() const;

        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        );

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
