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
    Foam::radiation::noReflection

Description
    Dummy scatter model for 'none'

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_noReflection_H
#define radiation_noReflection_H

#include "reflectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class noReflection Declaration
\*---------------------------------------------------------------------------*/

class noReflection
:
    public reflectionModel
{

public:

    //- Runtime type information
    TypeName("none");

    //- Construct from components
    noReflection(const dictionary& dict, const fvMesh& mesh);

    //- Destructor
    virtual ~noReflection() = default;


    // Member Functions

        //- Return reflectivity from medium1 to medium2 and a incident angle.
        virtual scalar rho
        (
            const scalar incidentAngle
        ) const;


        //- Return reflected direction
        virtual vector R(const vector& incident, const vector& n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
