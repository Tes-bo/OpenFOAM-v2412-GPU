/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::PlessisMasliyahDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle-drag model wherein drag forces (per unit carrier-fluid
    velocity) are dynamically computed based on the Du Plessis-Masliyah
    drag model.

    \f[
        \mathrm{F}_\mathrm{D} =
            \left(\mathrm{A}\, (1-\alpha_c) + \mathrm{B}\, \mathrm{Re}\right)
            \frac{(1-\alpha_c)\, \mu_c}{\alpha_c^2\, d_p^2}
    \f]
    with

    \f[
        A = \frac{26.8\, \alpha_c^2}
        {
            \alpha_p^{2/3}
            (1 - \alpha_p^{1/3})
            (1 - \alpha_p^{2/3})
        }
    \f]

    \f[
        \mathrm{B} = \frac{\alpha_c^2}{\left( 1 - \alpha_p^{2/3} \right)^2}
    \f]

    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    where
    \vartable
        \mathrm{F}_\mathrm{D} | Drag force per carrier-fluid velocity [kg/s]
        \mathrm{Re}_p         | Particle Reynolds number
        \mu_c    | Dynamic viscosity of carrier at the cell occupying particle
        d_p                   | Particle diameter
        \rho_c   | Density of carrier at the cell occupying particle
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
        \alpha_c              | Volume fraction of carrier fluid
        \alpha_p              | Volume fraction of particles
    \endvartable

    References:
    \verbatim
        Standard model (tag:P), (Eq. 34-36):
            Du Plessis, J. P. (1994).
            Analytical quantification of coefficients in the
            Ergun equation for fluid friction in a packed bed.
            Transport in porous media, 16(2), 189-207.
            DOI:10.1007/BF00617551
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        particleForces
        {
            PlessisMasliyahDrag
            {
                alphac    <alphacName>; // e.g. alpha.air
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type | Reqd | Deflt
      type      | Type name: PlessisMasliyahDrag    | word | yes  | -
      alphac    | Name of carrier fluid             | word | yes  | -
    \endtable

Note
  - \f$\mathrm{F}_\mathrm{D}\f$ is weighted with the particle mass/density
    at the stage of a function return, so that it can later be normalised
    with the effective mass, if necessary (e.g. when using virtual-mass forces).

SourceFiles
    PlessisMasliyahDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef PlessisMasliyahDragForce_H
#define PlessisMasliyahDragForce_H

#include "ParticleForce.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PlessisMasliyahDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PlessisMasliyahDragForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Reference to the carrier volume fraction field
        const volScalarField& alphac_;


public:

    //- Runtime type information
    TypeName("PlessisMasliyahDrag");


    // Constructors

        //- Construct from mesh
        PlessisMasliyahDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        PlessisMasliyahDragForce(const PlessisMasliyahDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new PlessisMasliyahDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const PlessisMasliyahDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~PlessisMasliyahDragForce() = default;


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
    #include "PlessisMasliyahDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
