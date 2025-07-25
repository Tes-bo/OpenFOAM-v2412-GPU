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
    Foam::Bezier

Description
    Calculation of adjoint based sensitivities for Bezier control points

SourceFiles
    Bezier.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Bezier_H
#define Foam_Bezier_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "volPointInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "deltaBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Bezier Declaration
\*---------------------------------------------------------------------------*/

class Bezier
{
private:

    // Private Member Functions

        //- No copy construct
        Bezier(const Bezier&) = delete;

        //- No copy assignment
        void operator=(const Bezier&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;
        const dictionary& dict_;

        label nBezier_;
        PtrList<pointTensorField> dxidXj_;

        boolList confineXmovement_;
        boolList confineYmovement_;
        boolList confineZmovement_;
        boolListList confineMovement_;
        labelList activeDesignVariables_;


public:

    //- Runtime type information
    TypeName("Bezier");


    // Constructors

        //- Construct from components
        Bezier(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~Bezier() = default;


    // Member Functions

        //- Number of Bezier control points
        label nBezier() const;

        //- dx/db tensor for all control points
        PtrList<pointTensorField>& dxidXj();

        //- Confine x movement
        const boolList& confineXmovement() const;

        //- Confine y movement
        const boolList& confineYmovement() const;

        //- Confine z movement
        const boolList& confineZmovement() const;

        //- Info about confining movement in all directions
        const boolListList& confineMovement() const;

        //- Compute derivative of the normal vector for a Bezier
        //- parameterized patch
        tmp<tensorField> dndbBasedSensitivities
        (
            const label patchI,
            const label cpI,
            bool returnDimensionedNormalSens = true
        ) const;

        //- Compute derivative of the normal vector for a Bezier parameterized
        //- patch wrt a given component (direction) of the Bezier control point
        tmp<vectorField> dndbBasedSensitivities
        (
            const label patchI,
            const label cpI,
            const label idir,
            bool returnDimensionedNormalSens = true
        ) const;

        //- dxdb tensor for a Bezier parameterized patch
        tmp<tensorField> dxdbFace
        (
            const label patchI,
            const label cpI,
            bool useChainRule = true
        ) const;

        //- dxdb vector for a Bezier parameterized patch and a specific control point component
        tmp<vectorField> dxdbFace
        (
            const label patchI,
            const label cpI,
            const label idir,
            bool useChainRule = true
        ) const;

        //- For a given (global) face ID, return the change of the face points
        tensorField facePoints_d
        (
            const label globalFaceI,
            const label cpI
        ) const;

        //- For a given (global) face ID, return the change of the face points
        //- if a given component of the CP changes
        vectorField facePoints_d
        (
            const label globalFaceI,
            const label cpI,
            const label idir
        ) const;

        //- Return active design variables
        const labelList& getActiveDesignVariables() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
