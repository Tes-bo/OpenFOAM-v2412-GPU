/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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
    Foam::DistortedSphereDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle-drag model wherein drag forces (per unit carrier-fluid velocity)
    are dynamically computed by using \c sphereDrag model; however, are
    corrected for particle distortion by linearly varying the drag between of
    a sphere (i.e. \c sphereDrag) and a value of 1.54 corresponding to a disk.

    \f[
        \mathrm{F}_\mathrm{D} =
            \frac{3}{4}
            \frac{\mu_c\,\mathrm{C}_\mathrm{D}\,\mathrm{Re}_p}{\rho_p \, d_p^2}
    \f]
    with

    \f[
        \mathrm{C}_\mathrm{D} =
            \mathrm{C}_{\mathrm{D, sphere}} \left( 1 + 2.632 y \right)
    \f]

    where
    \vartable
        \mathrm{F}_\mathrm{D}  | Drag force per carrier-fluid velocity [kg/s]
        \mathrm{C}_\mathrm{D}  | Particle drag coefficient
        \mathrm{C}_{\mathrm{D, sphere}} | Sphere drag coefficient
        \mathrm{Re}_p          | Particle Reynolds number
        \rho_p                 | Particle mass density
        d_p                    | Particle diameter
        y    | Level of distortion determined by other models internally
    \endvartable

    Constraints:
    - Applicable to particles with a spatially homogeneous distribution.
    - \f$ 1 \geq y \geq 0 \f$

    References:
    \verbatim
        Standard model:
            Putnam, A. (1961).
            Integratable form of droplet drag coefficient.
            ARS Journal, 31(10), 1467-1468.

        Standard model (tag:AOB):
            Amsden, A. A., O'Rourke, P. J., & Butler, T. D. (1989).
            KIVA-II: A computer program for chemically
            reactive flows with sprays (No. LA-11560-MS).
            Los Alamos National Lab.(LANL), Los Alamos, NM (United States).
            DOI:10.2172/6228444

        Expression correcting drag for particle distortion (tag:LMR):
            Liu, A. B., Mather, D., & Reitz, R. D. (1993).
            Modeling the effects of drop drag
            and breakup on fuel sprays.
            SAE Transactions, 83-95.
            DOI:10.4271/930072
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        particleForces
        {
            distortedSphereDrag;
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                       | Type | Reqd | Deflt
      type         | Type name: distortedSphereDrag    | word | yes  | -
    \endtable

Note
  - \f$\mathrm{F}_\mathrm{D}\f$ is weighted with the particle mass
    at the stage of a function return, so that it can later be normalised
    with the effective mass, if necessary (e.g. when using virtual-mass forces).

See also
  - Foam::SphereDragForce

SourceFiles
    DistortedSphereDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef DistortedSphereDragForce_H
#define DistortedSphereDragForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class DistortedSphereDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DistortedSphereDragForce
:
    public ParticleForce<CloudType>
{
    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;


public:

    //- Runtime type information
    TypeName("distortedSphereDrag");


    // Constructors

        //- Construct from mesh
        DistortedSphereDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        DistortedSphereDragForce(const DistortedSphereDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new DistortedSphereDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const DistortedSphereDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~DistortedSphereDragForce() = default;


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
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
    #include "DistortedSphereDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
