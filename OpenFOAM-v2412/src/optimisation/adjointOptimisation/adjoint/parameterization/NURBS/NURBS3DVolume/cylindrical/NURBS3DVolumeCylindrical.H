/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::NURBS3DVolumeCylindrical

Description
    Calculation of adjoint based sensitivities at wall faces

SourceFiles
    NURBS3DVolumeCylindrical.C

\*---------------------------------------------------------------------------*/

#ifndef NURBS3DVolumeCylindrical_H
#define NURBS3DVolumeCylindrical_H

#include "NURBS3DVolume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class NURBS3DVolumeCylindrical Declaration
\*---------------------------------------------------------------------------*/

class NURBS3DVolumeCylindrical
:
    public NURBS3DVolume
{
protected:

    // Protected Data

        //- Translation vector
        vector origin_;


    // Protected Member Functions

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
        NURBS3DVolumeCylindrical(const NURBS3DVolumeCylindrical&) = delete;

        //- No copy assignment
        void operator=(const NURBS3DVolumeCylindrical&) = delete;


public:

    //- Runtime type information
    TypeName("cylindrical");


    // Constructors

        //- Construct from components
        NURBS3DVolumeCylindrical
        (
            const dictionary& dict,
            const fvMesh& mesh,
            bool computeParamCoors = true
        );


    //- Destructor
    virtual ~NURBS3DVolumeCylindrical() = default;


    // Access

        //- Get the origin vector
        inline const vector& getOrigin() const
        {
            return origin_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
