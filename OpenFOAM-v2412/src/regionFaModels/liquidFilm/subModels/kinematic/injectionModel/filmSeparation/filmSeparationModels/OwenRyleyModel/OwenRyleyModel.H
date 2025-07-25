/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::filmSeparationModels::OwenRyleyModel

Description
    Computes film-separation properties from round edges for full separation
    (Owen & Ryley, 1983). The model assesses the film curvature based on the
    mesh geometry and flow field, and then calculates a force balance of the
    following form:

    \f[
        F_{net} = F_{inertial} + F_{body} + F_{surface}
    \f]

    with:

    \f[
        F_{inertial} = -\frac{4}{3} \frac{h \, \rho \, |\mathbf{u}|^2}{R_1}
    \f]

    \f[
        F_{body} =
            -\frac{\rho |\mathbf{g}| (R_1^2 - R_2^2) cos(\theta)}{2 \, R_1}
    \f]

    \f[
        F_{surface} = \frac{\sigma}{R_2}
    \f]

    where:

    \vartable
      F_{net}       | Magnitude of net force on the film
      F_{inertial}  | Magnitude of inertial forces on the film
      F_{body}      | Magnitude of body forces on the film
      F_{surface}   | Magnitude of surface forces on the film
      h             | Film thickness
      \rho          | Film density
      \sigma        | Film surface tension
      \mathbf{u}    | Film velocity
      \mathbf{g}    | Gravitational accelaration
      \theta | Approximate angle between gravity vector and outflow direction
      R_1           | Radius of curvature
      R_2           | Sum of the radius of curvature and film thickness
    \endvartable

    The film is considered separated when \f$F_{net}<0\f$; otherwise, it
    remains attached.

    Reference:
    \verbatim
    Governing equations (tag:OR):
        Owen, I., & Ryley, D. J. (1983).
        The flow of thin liquid films around corners.
        International journal of multiphase flow, 11(1), 51-62.
        DOI:10.1016/0301-9322(85)90005-9
    \endverbatim

Usage
    Minimal example in boundary-condition files:
    \verbatim
    filmSeparationCoeffs
    {
        // Mandatory entries
        model              OwenRyley;

        // Optional entries
        fThreshold         <scalar>;
        definedPatchRadii  <scalar>;
        deltaByR1Min       <scalar>;
        minInvR1           <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      model        | Model name: OwenRyley              | word | yes  | -
      fThreshold   | Threshold force for separation     | scalar | no | 1e-8
      definedPatchRadii | Patch radius                  | scalar | no | 0
      deltaByR1Min | Minimum gravity driven film thickness | scalar | no | 0
      minInvR1     | Minimum inv R1 for separation      | scalar | no | 5
    \endtable

Note
  - The assumptions underlying the derivation of the model indicate that
    the model is better suited for round corners than for sharp corners.

SourceFiles
    OwenRyleyModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_filmSeparationModels_OwenRyleyModel_H
#define Foam_filmSeparationModels_OwenRyleyModel_H

#include "filmSeparationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace filmSeparationModels
{

/*---------------------------------------------------------------------------*\
                    Class OwenRyleyModel Declaration
\*---------------------------------------------------------------------------*/

class OwenRyleyModel
:
    public filmSeparationModel
{
    // Private Data

        //- Threshold force for separation
        scalar fThreshold_;

        //- Patch radius
        scalar definedPatchRadii_;

        //- Minimum gravity driven film thickness (non-dimensionalised h/R1)
        scalar minHbyR1_;

        //- Minimum inv R1 of film velocity for film separation
        scalar minInvR1_;

        //- Gradient of surface normals
        areaTensorField gradNHat_;


    // Private Member Functions

        //- Calculate inverse of (local) radius of curvature of film velocity
        tmp<areaScalarField> calcInvR1(const areaVectorField& U) const;

        //- Calculate the cosine of the angle between gravity vector and
        //- cell-out flow direction
        tmp<scalarField> calcCosAngle
        (
            const edgeScalarField& phi,
            const scalarField& invR1
        ) const;

        //- Calculate the magnitude of net force on the film
        tmp<scalarField> netForce() const;


public:

    //- Runtime type information
    TypeName("OwenRyley");


    // Constructors

        //- Construct from the base film model and dictionary
        OwenRyleyModel
        (
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Destructor
    virtual ~OwenRyleyModel() = default;


    // Member Functions

    // Access

        // Return access to minimum dimensionless gravity driven film thickness
        scalar minHbyR1() const noexcept { return minHbyR1_; }

        // Return access to minimum inv R1 for film separation
        scalar minInvR1() const noexcept { return minInvR1_; }


    // Evaluation

        //- Calculate the mass ratio of film separation
        virtual tmp<scalarField> separatedMassRatio() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmSeparationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
