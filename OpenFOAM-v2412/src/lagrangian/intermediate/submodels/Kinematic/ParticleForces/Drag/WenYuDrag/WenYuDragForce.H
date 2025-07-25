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
    Foam::WenYuDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle-drag model wherein drag forces (per unit carrier-fluid
    velocity) are dynamically computed based on the Wen-Yu drag model.

    \f[
        \mathrm{F}_\mathrm{D} =
            \frac{3}{4}
            \frac{(1 - \alpha_c) \, \mu_c \, \alpha_c \, \mathrm{Re}_p }{d_p^2}
            \mathrm{C}_\mathrm{D} \, \alpha_c^{-2.65}
    \f]
    with

    \f[
        \mathrm{C}_\mathrm{D} =
            \frac{24}{\alpha_c \, \mathrm{Re}_p}
            \left(1 + \frac{1}{6}(\alpha_c \, \mathrm{Re}_p)^{2/3} \right)
            \quad \mathrm{if} \quad \alpha_c \, \mathrm{Re}_p < 1000
    \f]
    \f[
        \mathrm{C}_\mathrm{D} =
            0.44 \quad \mathrm{if} \quad \alpha_c \, \mathrm{Re}_p \geq 1000
    \f]
    and
    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    where
    \vartable
        \mathrm{F}_\mathrm{D} | Drag force per carrier-fluid velocity [kg/s]
        \mathrm{C}_\mathrm{D} | Particle drag coefficient
        \mathrm{Re}_p         | Particle Reynolds number
        \mu_c    | Dynamic viscosity of carrier at the cell occupying particle
        d_p                   | Particle diameter
        \rho_c   | Density of carrier at the cell occupying particle
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
        \alpha_c              | Volume fraction of the carrier fluid
    \endvartable

    References:
    \verbatim
        Standard model:
            Wen, C. Y., & Yu, Y. H., (1966).
            Mechanics of fluidization.
            Chem. Eng. Prog. Symp. Ser. 62, 100-111.

        Drag-coefficient model:
            Schiller, L., & Naumann, A. (1935).
            Über die grundlegenden Berechnungen bei der Schwerkraftaufbereitung.
            Z. Ver. Dtsch. Ing., 77: 318–326.

        Expressions (tags:ZZB, GLSLR), (Eq.13-14, Table 3):
            Zhou, L., Zhang, L., Bai, L., Shi, W.,
            Li, W., Wang, C., & Agarwal, R. (2017).
            Experimental study and transient CFD/DEM simulation in
            a fluidized bed based on different drag models.
            RSC advances, 7(21), 12764-12774.
            DOI:10.1039/C6RA28615A

            Gao, X., Li, T., Sarkar, A., Lu, L., & Rogers, W. A. (2018).
            Development and validation of an enhanced filtered drag model
            for simulating gas-solid fluidization of Geldart A particles
            in all flow regimes.
            Chemical Engineering Science, 184, 33-51.
            DOI:10.1016/j.ces.2018.03.038
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        particleForces
        {
            WenYuDrag
            {
                alphac    <alphacName>; // e.g. alpha.air
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                       | Type | Reqd | Deflt
      type         | Type name: WenYuDrag              | word | yes  | -
      alphac       | Name of carrier fluid             | word | yes  | -
    \endtable

Note
  - \f$\mathrm{F}_\mathrm{D}\f$ is weighted with the particle mass/density
    at the stage of a function return, so that it can later be normalised
    with the effective mass, if necessary (e.g. when using virtual-mass forces).

SourceFiles
    WenYuDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef WenYuDragForce_H
#define WenYuDragForce_H

#include "ParticleForce.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class WenYuDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WenYuDragForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Reference to the carrier volume fraction field
        const volScalarField& alphac_;


    // Private Member Functions

        //- Drag coefficient multiplied by volume fraction and Reynolds number
        scalar CdRe(const scalar alphacRe) const;


public:

    //- Runtime type information
    TypeName("WenYuDrag");


    // Constructors

        //- Construct from mesh
        WenYuDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        WenYuDragForce(const WenYuDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new WenYuDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const WenYuDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~WenYuDragForce() = default;


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
    #include "WenYuDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
