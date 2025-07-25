/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::seriesProfile

Description
    Profile model where polar lift and drag coefficients are computed
    as sum of trigonometric series by using an input angle of attack.

    \f[
        C_d = \sum_i\left(C_{dCoeff} \cos(i \alpha) \right)
    \f]

    \f[
        C_l = \sum_i\left(C_{lCoeff} \sin(i \alpha) \right)
    \f]

    where
    \vartable
      C_d    | Drag coefficient
      C_l    | Lift coefficient
      \alpha | Angle of attack      [deg]
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSource1
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (runtime modifiable)
        profiles
        {
            // Mandatory entries (runtime modifiable)
            <profile1>
            {
                // Mandatory entries (runtime modifiable)
                type     series;
                CdCoeffs (coeff1 coeff2 ... coeffN);
                ClCoeffs (coeff1 coeff2 ... coeffN);
            }

            // Mandatory/Optional (inherited) entries
            ...
        }
    }

    where the entries mean:
    \table
      Property   | Description                   | Type       | Reqd | Dflt
      CdCoeffs   | List of drag coefficients     | scalarList | yes  | -
      ClCoeffs   | List of lift coefficients     | scalarList | yes  | -
    \endtable

Note
  - Angle of attack is internally converted from [deg] to [rad].

See also
  - Foam::fv::rotorDiskSource
  - Foam::profileModel
  - Foam::seriesProfile

SourceFiles
    seriesProfile.C

\*---------------------------------------------------------------------------*/

#ifndef seriesProfile_H
#define seriesProfile_H

#include "profileModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class seriesProfile Declaration
\*---------------------------------------------------------------------------*/

class seriesProfile
:
    public profileModel
{
protected:

    // Protected Data

        //- List of drag coefficient values
        List<scalar> CdCoeffs_;

        //- List of lift coefficient values
        List<scalar> ClCoeffs_;


    // Protected Member Functions

        // Evaluation

            //- Drag
            scalar evaluateDrag
            (
                const scalar& xIn,
                const List<scalar>& values
            ) const;

            //- Lift
            scalar evaluateLift
            (
                const scalar& xIn,
                const List<scalar>& values
            ) const;


public:

    //- Runtime type information
    TypeName("series");


    // Constructors

        //- Constructor from dictionary and model name
        seriesProfile(const dictionary& dict, const word& modelName);

        //- No copy construct
        seriesProfile(const seriesProfile&) = delete;

        //- No copy assignment
        void operator=(const seriesProfile&) = delete;


    //- Destructor
    ~seriesProfile() = default;


    // Member Functions

        //- Return the Cd and Cl for a given angle-of-attack
        virtual void Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
