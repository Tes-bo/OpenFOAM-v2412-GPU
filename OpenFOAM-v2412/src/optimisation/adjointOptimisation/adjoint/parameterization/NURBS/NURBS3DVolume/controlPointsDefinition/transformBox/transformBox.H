/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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
    Foam::transformBox

Description
    Defines control points based on given lower and upper boundaries for the
    box and additional translation, rotation and scaling vectors.
    The latter can be copied over from the Paraview clip or transform filters,
    to facilitate the definition of a more general box.

SourceFiles
    transformBox.C

\*---------------------------------------------------------------------------*/

#ifndef transformBox_H
#define transformBox_H

#include "controlPointsDefinition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class transformBox Declaration
\*---------------------------------------------------------------------------*/

class transformBox
:
    public controlPointsDefinition
{

protected:

    // Protected Member Functions

        //- Compute control points
        virtual void computeControlPoints();


private:

    // Private Member Functions

        //- No copy construct
        transformBox(const transformBox&) = delete;

        //- No copy assignment
        void operator=(const transformBox&) = delete;


public:

    //- Runtime type information
    TypeName("transformBox");


    // Constructors

        //- Construct from components
        transformBox(NURBS3DVolume& box);


    //- Destructor
    virtual ~transformBox() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
