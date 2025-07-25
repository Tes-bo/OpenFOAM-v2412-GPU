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
    Foam::lookupProfile

Description
    Profile model where polar lift and drag coefficients are linearly
    interpolated from a polar table by using an input angle of attack.

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
                type     lookup;
                data
                (
                    (AOA1 Cd1 Cl2)
                    (AOA2 Cd2 Cl2)
                    ([0]  [1] [2])
                    ...
                    (AOAN CdN CdN)
                );
            }

            // Mandatory/Optional (inherited) entries
            ...
        }
    }

    where the entries mean:
    \table
      Property   | Description              | Type   | Reqd | Dflt
      [0] AOA    | Angle of attack          | scalar | yes  | -
      [1] Cd     | Drag coefficient corresponding to the angle attack <!--
                 -->                        | scalar | yes  | -
      [2] Cl     | Lift coefficient corresponding to the angle attack <!--
                 -->                        | scalar | yes  | -
    \endtable

Note
  - Angle of attack is internally converted from [deg] to [rad].

See also
  - Foam::fv::rotorDiskSource
  - Foam::profileModel
  - Foam::seriesProfile

SourceFiles
    lookupProfile.C

\*---------------------------------------------------------------------------*/

#ifndef lookupProfile_H
#define lookupProfile_H

#include "profileModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class lookupProfile Declaration
\*---------------------------------------------------------------------------*/

class lookupProfile
:
    public profileModel
{
protected:

    // Protected Data

        //- List of angle-of-attack values [deg] on input, converted to [rad]
        List<scalar> AOA_;

        //- List of drag coefficient values
        List<scalar> Cd_;

        //- List of lift coefficient values
        List<scalar> Cl_;


    // Protected Member Functions

        //- Return the interpolation indices and gradient
        void interpolateWeights
        (
            const scalar& xIn,
            const List<scalar>& values,
            label& i1,
            label& i2,
            scalar& ddx
        ) const;


public:

    //- Runtime type information
    TypeName("lookup");


    // Constructors

        //- Constructor from dictionary and model name
        lookupProfile(const dictionary& dict, const word& modelName);

        //- No copy construct
        lookupProfile(const lookupProfile&) = delete;

        //- No copy assignment
        void operator=(const lookupProfile&) = delete;


    //- Destructor
    ~lookupProfile() = default;


    // Member Functions

        //- Return the Cd and Cl for a given angle-of-attack
        virtual void Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
