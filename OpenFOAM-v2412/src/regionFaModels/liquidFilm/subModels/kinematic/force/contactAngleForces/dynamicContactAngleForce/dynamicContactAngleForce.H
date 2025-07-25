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
    Foam::regionModels::areaSurfaceFilmModels::dynamicContactAngleForce

Description
    Film-speed or film-temperature dependent
    contact-angle force with a stochastic perturbation.

    The contact angle in degrees is specified as a \c Foam::Function1 type,
    to enable the use of, e.g. \c constant, \c polynomial, \c table values
    and the stochastic perturbation obtained from a
    \c Foam::distributionModels::distributionModel.

Usage
    Minimal example:
    \verbatim
    forces
    (
        dynamicContactAngle
    );

    dynamicContactAngleForceCoeffs
    {
        // Mandatory entries
        distribution    <subDict>;

        // Conditional entries

            // Option-1
            Utheta           <Function1<scalar>>;

            // Option-2
            Ttheta           <Function1<scalar>>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property              | Description              | Type  | Reqd | Deflt
      dynamicContactAngle   | Type name                | word  | yes  | -
      Utheta      | Contact angle as a function of film speed <!--
            -->                        | \<Function1\<scalar\> | choice | -
      Ttheta      | Contact angle as a function of film temperature <!--
            -->                        | \<Function1\<scalar\> | choice | -
      distribution  | Probability distribution model  | subDict | yes | -
    \endtable

    The inherited entries are elaborated in:
      - \link contactAngleForce.H \endlink
      - \link Function1.H \endlink
      - \link distributionModel.H \endlink

SourceFiles
    dynamicContactAngleForce.C

\*---------------------------------------------------------------------------*/

#ifndef areaSurfaceFilmModels_dynamicContactAngleForce_H
#define areaSurfaceFilmModels_dynamicContactAngleForce_H

#include "contactAngleForce.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class Function1;
class distributionModel;

namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                  Class dynamicContactAngleForce Declaration
\*---------------------------------------------------------------------------*/

class dynamicContactAngleForce
:
    public contactAngleForce
{
    // Private Data

        //- Contact angle as a function of film speed
        autoPtr<Function1<scalar>> U_vs_thetaPtr_;

        //- Contact angle as a function of film temperature
        autoPtr<Function1<scalar>> T_vs_thetaPtr_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<distributionModel> distribution_;


    // Private Member Functions

        //- No copy construct
        dynamicContactAngleForce(const dynamicContactAngleForce&) = delete;

        //- No copy assignment
        void operator=(const dynamicContactAngleForce&) = delete;


protected:

        //- Return the contact angle field
        virtual tmp<areaScalarField> theta() const;


public:

    //- Runtime type information
    TypeName("dynamicContactAngle");


    // Constructors

        //- Construct from surface film model
        dynamicContactAngleForce
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dynamicContactAngleForce() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
