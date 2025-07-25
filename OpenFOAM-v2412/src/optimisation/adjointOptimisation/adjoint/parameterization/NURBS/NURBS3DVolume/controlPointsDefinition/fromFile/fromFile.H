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
    Foam::fromFile

Description
    Reads control points from constant/controlPoints

SourceFiles
    fromFile.C

\*---------------------------------------------------------------------------*/

#ifndef fromFile_H
#define fromFile_H

#include "controlPointsDefinition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fromFile Declaration
\*---------------------------------------------------------------------------*/

class fromFile
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
        fromFile(const fromFile&) = delete;

        //- No copy assignment
        void operator=(const fromFile&) = delete;


public:

    //- Runtime type information
    TypeName("fromFile");


    // Constructors

        //- Construct from components
        fromFile(NURBS3DVolume& box);


    //- Destructor
    virtual ~fromFile() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
