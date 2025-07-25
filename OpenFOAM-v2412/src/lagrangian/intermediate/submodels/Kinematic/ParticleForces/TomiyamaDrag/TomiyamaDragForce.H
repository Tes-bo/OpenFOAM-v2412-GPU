/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::TomiyamaDragForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle-drag model wherein drag forces (per unit carrier-fluid
    velocity) are dynamically computed using empirical expressions based on
    their level of contamination.

    \f[
        \mathrm{F}_\mathrm{D} =
            \frac{3}{4}
            \frac{\mu_c\,\mathrm{C}_\mathrm{D}\,\mathrm{Re}_p}{\rho_p \, d_p^2}
    \f]

    For pure systems:
    \f[
        \mathrm{C}_\mathrm{D} =
            \max\left[
                \min\left(\frac{16}{Re}(1+0.15Re^{0.687}), \frac{48}{Re}\right),
                \frac{8}{3}\frac{Eo}{Eo+4}
            \right]
    \f]

    For slightly contaminated systems:
    \f[
        \mathrm{C}_\mathrm{D} =
            \max\left[
                \min\left(\frac{24}{Re}(1+0.15Re^{0.687}), \frac{72}{Re}\right),
                \frac{8}{3}\frac{Eo}{Eo+4}
            \right]
    \f]

    For fully contaminated systems:
    \f[
        \mathrm{C}_\mathrm{D} =
            \max\left[
                \frac{24}{Re}(1+0.15Re^{0.687}),
                \frac{8}{3}\frac{Eo}{Eo+4}
            \right]
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
        Eo       | Eotvos number
    \endvartable

    Constraints:
    - Applicable to bubbles with a spatially homogeneous distribution.

    References:
    \verbatim
        Tomiyama, A., Kataoka, I., Zun, I., & Sakaguchi, T. (1998).
        Drag coefficients of single bubbles under normal and micro gravity
        conditions.
        JSME International Journal Series B
        Fluids and Thermal Engineering, 41(2), 472-479.
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        particleForces
        {
            tomiyamaDrag
            {
                // Mandatory entries
                sigma         <scalar>;
                contamination <word>;  // pure | slight | full
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                       | Type   | Reqd | Deflt
      type          | Type name: tomiyamaDrag           | word   | yes  | -
      sigma         | Surface tension                   | scalar | yes  | -
      contamination | Contamination type                | word   | yes  | -
    \endtable

    Options for the \c contamination entry:
    \verbatim
      pure          | Pure systems
      slight        | Slightly contaminated systems
      full          | Fully contaminated systems
    \endverbatim

Note
  - \f$\mathrm{F}_\mathrm{D}\f$ is weighted with the particle mass/density
    at the stage of a function return, so that it can later be normalised
    with the effective mass, if necessary (e.g. when using virtual-mass forces).

SourceFiles
    TomiyamaDragForce.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaDragForce_H
#define TomiyamaDragForce_H

#include "ParticleForce.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class TomiyamaDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TomiyamaDragForce
:
    public ParticleForce<CloudType>
{
public:

    // Public Enumerations

        //- Options for the contamination types
        enum contaminationType : char
        {
            PURE = 0,  //!< "Pure systems"
            SLIGHT,    //!< "Slightly contaminated systems"
            FULL       //!< "Fully contaminated systems"
        };

        //- Names for the contaminationType options
        static const Enum<contaminationType> contaminationTypeNames;


private:

    // Private Data

        //- Surface tension
        const scalar sigma_;

        //- Contamination type option
        const contaminationType contaminationType_;


    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar CdRe(const scalar Re) const;


public:

    //- Runtime type information
    TypeName("TomiyamaDrag");


    // Constructors

        //- Construct from mesh
        TomiyamaDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy construct
        TomiyamaDragForce(const TomiyamaDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new TomiyamaDragForce<CloudType>(*this)
            );
        }

        //- No copy assignment
        void operator=(const TomiyamaDragForce<CloudType>&) = delete;


    //- Destructor
    virtual ~TomiyamaDragForce() = default;


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
    #include "TomiyamaDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
