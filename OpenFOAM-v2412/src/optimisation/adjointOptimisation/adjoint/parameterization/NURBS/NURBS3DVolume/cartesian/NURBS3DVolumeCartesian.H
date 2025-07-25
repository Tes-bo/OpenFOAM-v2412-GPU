/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::NURBS3DVolumeCartesian

Description
    Calculation of adjoint based sensitivities at wall faces

SourceFiles
    NURBS3DVolumeCartesian.C

\*---------------------------------------------------------------------------*/

#ifndef NURBS3DVolumeCartesian_H
#define NURBS3DVolumeCartesian_H

#include "NURBS3DVolume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class NURBS3DVolumeCartesian Declaration
\*---------------------------------------------------------------------------*/

class NURBS3DVolumeCartesian
:
    public NURBS3DVolume
{
protected:

    // Protected data

        //- Transform a point from its coordinate system to a cartesian system
        vector transformPointToCartesian(const vector& localCoordinates) const;

        //- Transformation tensor for dxdb, from local coordinate system to
        //- cartesian
        tensor transformationTensorDxDb(label globalPointIndex);

        //- Update coordinates in the local system based on the cartesian points
        void updateLocalCoordinateSystem(const vectorField& cartesianPoints);


private:

    // Private Member Functions

        //- No copy construct
        NURBS3DVolumeCartesian(const NURBS3DVolumeCartesian&) = delete;

        //- No copy assignment
        void operator=(const NURBS3DVolumeCartesian&) = delete;


public:

    //- Runtime type information
    TypeName("cartesian");


    // Constructors

        //- Construct from components
        NURBS3DVolumeCartesian
        (
            const dictionary& dict,
            const fvMesh& mesh,
            bool computeParamCoors = true
        );


    //- Destructor
    virtual ~NURBS3DVolumeCartesian() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
