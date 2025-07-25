/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::contactAngleForce

Description
    Base-class for film contact angle force models.

    The effect of the contact angle force can be ignored over a specified
    distance from patches.

    The contact-angle force is implemented as follows:
    \f[
        \vec{f}_\theta =
            \beta \frac{\sigma (1 - cos(\theta))}{\Delta_{cl}} \vec{n}_{cl}
            \frac{1}{\rho_f|\vec{S}_f|}
          - \vec{f}_{cl}
    \f]

    with if \f$ h > h_{crit} \f$:

    \f[
        \vec{f}_{cl} = \nu_f \frac{\vec{U}_f}{h_{crit}}
    \f]

    where
    \vartable
        \vec{f}_\theta | Contact-angle force            [N/m^2/rho]
        \beta          | Empirical constant             [-]
        \sigma         | Liquid-air interfacial surface tension [kg/s^2]
        \theta         | Contact angle                  [rad]
        \Delta_{cl}    | Inverse width in direction normal to contact line [1/m]
        \vec{n}_{cl}   | Direction normal to contact line [-]
        \rho_f         | Film density                   [kg/m^3]
        \vec{S}_f      | Face-area vector               [m^2]
        \vec{f}_{cl}   | Contact-line movement force    [N/m^2/rho]
        \nu_f          | Film kinematic viscosity       [m^2/s]
        \vec{U}_f      | Film velocity                  [m/s]
        h_{crit}       | Critical film height           [m]
    \endvartable

    Reference:
    \verbatim
        Governing equations (tag:MHDX):
            Meredith, K. V., Heather, A., De Vries, J., & Xin, Y. (2011).
            A numerical model for partially-wetted flow of thin liquid films.
            Computational Methods in Multiphase Flow VI, 70, 239.

        Contact line movement (tag:NDPC):
            Novák, M., Devaradja, R., Papper, J., & Černý, M. (2020).
            Efficient CFD methods for assessment of water management.
            In 20. Internationales Stuttgarter Symposium (pp. 151-170).
            Springer Vieweg, Wiesbaden.
    \endverbatim

Usage
    Minimal example:
    \verbatim
    {
        // Mandatory entries
        Ccf             <scalar>;

        // Optional entries
        hCrit           <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property         | Description               | Type  | Reqd | Deflt
      Ccf              | Empirical coefficient    | scalar | yes  | -
      hCrit            | Critical film height [m] | scalar | no   | GREAT
    \endtable

    The inherited entries are elaborated in:
      - \link force.H \endlink

SourceFiles
    contactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef areaSurfaceFilmModels_contactAngleForce_H
#define areaSurfaceFilmModels_contactAngleForce_H

#include "force.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class contactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class contactAngleForce
:
    public force
{
    // Private Data

        //- Coefficient applied to the contact angle force
        scalar Ccf_;

        //- Critical film thickness [m]
        scalar hCrit_;


    // Private Member Functions

        //- No copy construct
        contactAngleForce(const contactAngleForce&) = delete;

        //- No copy assignment
        void operator=(const contactAngleForce&) = delete;


protected:

        //- Return the contact angle field
        virtual tmp<areaScalarField> theta() const = 0;


public:

    //- Runtime type information
    TypeName("contactAngle");


    // Constructors

        //- Construct from surface film model
        contactAngleForce
        (
            const word& typeName,
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~contactAngleForce() = default;


    // Member Functions

        //- Correct
        virtual tmp<faVectorMatrix> correct(areaVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
