/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::NonSphereDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle-drag model wherein drag forces (per unit carrier-fluid
    velocity) are dynamically computed based on empirical expressions using
    a four-parameter general drag correlation for non-spherical particles.

    \f[
        \mathrm{F}_\mathrm{D} =
            \frac{3}{4}
            \frac{\mu_c\,\mathrm{C}_\mathrm{D}\,\mathrm{Re}_p}{\rho_p \, d_p^2}
    \f]
    with

    \f[
        \mathrm{C}_\mathrm{D} =
            \frac{24}{\mathrm{Re}_p} \left( 1 + A \, \mathrm{Re}_p^B \right)
          + \frac{C \, \mathrm{Re}_p}{D + \mathrm{Re}_p}
    \f]
    where

    \f[
        A = \exp(2.3288 - 6.4581\phi + 2.4486 \phi^2)
    \f]

    \f[
        B = 0.0964 + 0.5565\phi
    \f]

    \f[
        C = \exp(4.9050 - 13.8944\phi + 18.4222\phi^2 - 10.2599 \phi^3)
    \f]

    \f[
        D = \exp(1.4681 + 12.2584\phi - 20.7322\phi^2 + 15.8855\phi^3)
    \f]

    \f[
        \phi = \frac{A_p}{A_a}
    \f]

    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    where
    \vartable
        \mathrm{F}_\mathrm{D} | Drag force per carrier-fluid velocity [kg/s]
        \mathrm{C}_\mathrm{D} | Particle drag coefficient
        \mathrm{Re}_p         | Particle Reynolds number
        \rho_p                | Particle mass density
        \mu_c    | Dynamic viscosity of carrier at the cell occupying particle
        d_p                   | Particle diameter
        \rho_c   | Density of carrier at the cell occupying particle
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
        A_p        | Surface area of sphere with the same volume as the particle
        A_a           | Actual surface area of the particle
        \phi          | Ratio of surface areas
    \endvartable

    Constraints:
    - Applicable to particles with a spatially homogeneous distribution.
    - \f$ 1 \geq \phi > 0 \f$

    References:
    \verbatim
        Standard model (tag:HL), (Eq. 4,10-11):
            Haider, A., & Levenspiel, O. (1989).
            Drag coefficient and terminal velocity of
            spherical and nonspherical particles.
            Powder technology, 58(1), 63-70.
            DOI:10.1016/0032-5910(89)80008-7
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        particleForces
        {
            nonSphereDrag
            {
                phi    <phi>;
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type | Reqd | Deflt
      type      | Type name: nonSphereDrag          | word | yes  | -
      phi       | Ratio of surface area of sphere having same     <!--
                --> volume as particle to actual surface area of  <!--
                --> particle                      | scalar | yes  | -
    \endtable

Note
  - The drag coefficient model in (HL:Eq. 11) is good to within
    2 to 4 \% of RMS values from the corresponding experiment.
  - (HL:Eq. 12) also give a simplified model with greater error
    compared to results from the experiment, but since \c phi is
    presumed constant, Eq. 12 offers little benefit.
  - \f$\mathrm{F}_\mathrm{D}\f$ is weighted with the particle mass
    at the stage of a function return, so that it can later be normalised
    with the effective mass, if necessary (e.g. when using virtual-mass forces).

SourceFiles
    NonSphereDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef NonSphereDragForce_H
#define NonSphereDragForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class NonSphereDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NonSphereDragForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected Data

        //- Ratio of surface of sphere having same volume as particle to
        //- actual surface area of particle (0 < phi <= 1)
        scalar phi_;


        // Model coefficients

            scalar a_;

            scalar b_;

            scalar c_;

            scalar d_;


    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;


public:

    //- Runtime type information
    TypeName("nonSphereDrag");


    // Constructors

        //- Construct from mesh
        NonSphereDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        NonSphereDragForce(const NonSphereDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new NonSphereDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const NonSphereDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~NonSphereDragForce() = default;


    // Member Functions

        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NonSphereDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
