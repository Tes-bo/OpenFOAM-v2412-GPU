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
    Foam::radiation::FresnelLaser

Description
    Modified Fresnel reflection model.

    \verbatim
    Implementation of real-time multiple reflection and Fresnel absorption
    of FresnelLaser beam in keyhole.
    J. Phys. D: Appl. Phys. 39 (2006) 5372-5378
    \endverbatim

SourceFiles
    FresnelLaser.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_FresnelLaser_H
#define radiation_FresnelLaser_H

#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                       Class FresnelLaser Declaration
\*---------------------------------------------------------------------------*/

class FresnelLaser
:
    public reflectionModel
{

    // Private data

        //- Model constant
        scalar epsilon_;


public:

    //- Runtime type information
    TypeName("FresnelLaser");

    //- Construct from components
    FresnelLaser(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~FresnelLaser() = default;


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
