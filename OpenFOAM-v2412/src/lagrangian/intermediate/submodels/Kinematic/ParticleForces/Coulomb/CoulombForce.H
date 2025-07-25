/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::CoulombForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle electric force model involving the Coulomb force calculation.

    \f[
        \vec{F}_\mathrm{E} = q \, \vec{E}
    \f]

    where

    \vartable
        \vec{F}_\mathrm{E}  | Coulomb force                 [kg m s^{-2}]
        q                   | Electric charge of particle   [A s]
        \vec{E}             | Electric field    [kg m^2 A^{-1} s^{-3} m^{-1}]
    \endvartable

    References:
    \verbatim
        Governing equations (tag:YSSD):
            Ye, Q., Steigleder, T., Scheibe, A., & Domnick, J. (2002).
            Numerical simulation of the electrostatic powder coating process
            with a corona spray gun. Journal of Electrostatics, 54(2), 189-205.
            DOI:10.1016/S0304-3886(01)00181-4
    \endverbatim

Usage
    Minimal example by using \c constant/cloudProperties:
    \verbatim
    subModels
    {
        solution
        {
            interpolationSchemes
            {
                <Ename>     <interpolationScheme>;
            }
        }

        particleForces
        {
            Coulomb
            {
                q       <Function1<scalar>>;
                E       <word>;
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                       | Type | Reqd | Deflt
      type     | Type name: Coulomb                | word | yes  | -
      q        | Electric charge of particles | \<Function1\<scalar\> | yes | -
      E        | Name of electric field            | word | no   | E
    \endtable

Note
  - Particle electric charge can be input as a function of particle diameter:
    \verbatim
        particleForces
        {
            Coulomb
            {
                q   table
                (
                    // d [m]           q [C = A s]
                    (1e-06             -5.0e-11)
                    (1e-05             -1.0e-10)
                );
                E   electricPotential:E;
            }
        }
    \endverbatim

SourceFiles
    CoulombForce.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_CoulombForce_H
#define Foam_CoulombForce_H

#include "ParticleForce.H"
#include "Function1.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class CoulombForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CoulombForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Particle electric charge
        autoPtr<Function1<scalar>> qPtr_;

        //- Name of electric field
        const word Ename_;

        //- Electric-field interpolator
        mutable std::unique_ptr<interpolation<vector>> EInterpPtr_;


    // Private Member Functions

        //- Return requested volVectorField from the object registry
        //- or read+register the field to the object registry
        volVectorField& getOrReadField(const word& fieldName) const;


public:

    //- Runtime type information
    TypeName("Coulomb");


    // Constructors

        //- Construct from mesh
        CoulombForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Copy construct
        CoulombForce(const CoulombForce& gf);

        //- No copy assignment
        void operator=(const CoulombForce<CloudType>&) = delete;

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new CoulombForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CoulombForce() = default;


    // Member Functions

        // Evaluation

        //- Cache fields
        virtual void cacheFields(const bool store);

        //- Calculate the non-coupled force
        virtual forceSuSp calcNonCoupled
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
    #include "CoulombForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
