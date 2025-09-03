/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::pointToPointPlanarInterpolation

Description
    Interpolates between two sets of unstructured points using 2D Delaunay
    triangulation. Used in e.g. timeVaryingMapped bcs.

SourceFiles
    pointToPointPlanarInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointToPointPlanarInterpolation_H
#define Foam_pointToPointPlanarInterpolation_H

#include "FixedList.H"
#include "instantList.H"
#include "Pair.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pointToPointPlanarInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointToPointPlanarInterpolation
{
    // Private Data

        //- Perturbation factor
        scalar perturb_;

        //- Use nearest point only (avoids triangulation, projection)
        bool nearestOnly_;

        //- Cartesian reference coordinate system
        coordSystem::cartesian referenceCS_;

        //- Number of source points (for checking)
        label nPoints_;

        //- Interpolation addressing to face centres of underlying patch
        List<FixedList<label, 3>> nearestVertex_;

        //- Interpolation factors to face centres of underlying patch
        List<FixedList<scalar, 3>> nearestVertexWeight_;


    // Private Member Functions

        //- Calculate a local coordinate system from set of points
        static coordSystem::cartesian calcCoordinateSystem(const pointField&);

        //- Calculate addressing and weights
        void calcWeights
        (
            const pointField& sourcePoints,
            const pointField& destPoints
        );


public:

    // Declare name of the class and its debug switch
    ClassName("pointToPointPlanarInterpolation");


    // Generated Methods

        //- Copy construct
        pointToPointPlanarInterpolation
        (
            const pointToPointPlanarInterpolation&
        ) = default;

        //- Copy assignment
        pointToPointPlanarInterpolation& operator=
        (
            const pointToPointPlanarInterpolation&
        ) = default;


    // Constructors

        //- Construct from 3D locations.
        //  Determines local coordinate system
        //  from sourcePoints and maps onto that. If nearestOnly skips any
        //  local coordinate system and triangulation and uses nearest vertex
        //  only
        pointToPointPlanarInterpolation
        (
            const pointField& sourcePoints,
            const pointField& destPoints,
            const scalar perturb,
            const bool nearestOnly = false
        );

        //- Construct from coordinate system and locations.
        pointToPointPlanarInterpolation
        (
            const coordinateSystem& referenceCS,
            const pointField& sourcePoints,
            const pointField& destPoints,
            const scalar perturb
        );

        //- Construct from components
        pointToPointPlanarInterpolation
        (
            const scalar perturb,
            const bool nearestOnly,
            const coordinateSystem& referenceCS,
            const label sourceSize,
            List<FixedList<label, 3>>&& nearestVertex,
            List<FixedList<scalar, 3>>&& nearestVertexWeight
        );

        //- Construct and return a clone
        autoPtr<pointToPointPlanarInterpolation> clone() const
        {
            return autoPtr<pointToPointPlanarInterpolation>::New(*this);
        }


    // Member Functions

        //- Perturbation factor (for triangulation)
        scalar perturb() const noexcept
        {
            return perturb_;
        }

        //- Whether to use nearest point only (avoids triangulation, projection)
        bool nearestOnly() const noexcept
        {
            return nearestOnly_;
        }

        //- Return the Cartesian reference coordinate system
        const coordSystem::cartesian& referenceCS() const noexcept
        {
            return referenceCS_;
        }

        //- Number of source points
        label sourceSize() const noexcept
        {
            return nPoints_;
        }

        //- Number of target points
        label targetSize() const noexcept
        {
            return nearestVertex_.size();
        }

        //- Interpolation addressing to face centres of underlying patch
        const List<FixedList<label, 3>>& nearestVertex() const noexcept
        {
            return nearestVertex_;
        }

        //- Interpolation factors to face centres of underlying patch
        const List<FixedList<scalar, 3>>& nearestVertexWeight() const noexcept
        {
            return nearestVertexWeight_;
        }

        //- Helper: extract words of times
        static wordList timeNames(const instantList& times);

        //- Interpolate from field on source points to dest points
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& sourceFld) const;

        //- Interpolate from field on source points to dest points
        template<class Type>
        tmp<Field<Type>> interpolate(const tmp<Field<Type>>& tsource) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointToPointPlanarInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
