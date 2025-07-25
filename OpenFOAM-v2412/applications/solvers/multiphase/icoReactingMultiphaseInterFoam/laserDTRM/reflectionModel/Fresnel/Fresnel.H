/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::radiation::Fresnel

Description
    General Fresnel reflection model bewtween a dialectric and an absorbing
    medium.

    Radiative heat transfer. Micheal Modest. 3dr Edition. Chapter 2.5

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_Fresnel_H
#define radiation_Fresnel_H

#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class Fresnel Declaration
\*---------------------------------------------------------------------------*/

class Fresnel
:
    public reflectionModel
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Complex index of refraction for medium 1
        Pair<scalar> nk1_;

        //- Complex index of refraction for medium 2
        Pair<scalar> nk2_;


public:

    //- Runtime type information
    TypeName("Fresnel");

    //- Construct from components
    Fresnel(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~Fresnel() = default;


    // Member Functions

        //- Return reflected direction
        virtual vector R(const vector& incident, const vector& n) const;

        //- Return reflectivity from medium1 to medium2 and a incident angle.
        //  nk1 = (n1 - i k1) from medium 1.
        virtual scalar rho(const scalar incidentAngle) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
