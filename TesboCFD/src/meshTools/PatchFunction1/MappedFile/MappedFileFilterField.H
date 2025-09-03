/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::PatchFunction1Types::FilterField

Description
    The FilterField helper class provides a multi-sweep median filter
    for a Field of data associated with a geometric point cloud.

    The points can be freestanding or the faceCentres (or points)
    of a meshedSurface, for example.

    Using an initial specified search radius, the nearest point
    neighbours are gathered and addressing/weights are built for them.
    This currently uses an area-weighted, linear RBF interpolator
    with provision for quadratic RBF interpolator etc.

    After the weights and addressing are established,
    the evaluate() method can be called to apply a median filter
    to data fields, with a specified number of sweeps.

Note
    When handling large search radii and/or an extensive number of
    filter sweeps, compiling with openmp can yield some speedup.

SourceFiles
    MappedFileFilterField.C
    MappedFileFilterFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PatchFunction1Types_MappedFileFilterField_H
#define Foam_PatchFunction1Types_MappedFileFilterField_H

#include "primitiveFields.H"
#include "pointField.H"
#include "Enum.H"
#include "className.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                         Class FilterField Declaration
\*---------------------------------------------------------------------------*/

class FilterField
{
public:

     // Data Types

        //- Basis type
        enum RBF_type
        {
            RBF_linear,
            RBF_quadratic
        };


private:

    // Private Data

        //- Names for RBF basis types
        static const Enum<RBF_type> RBF_typeNames_;

        //- Addressing
        List<labelList> addressing_;

        //- Weights
        List<scalarField> weights_;


    // Private Member Functions

        //- Construct neighbour point weights and addressing based
        //- on specified search radius
        //
        // Search with local sphere
        //
        // RadiusSqrOp (index) -> radiusSqr
        // BasisFunctionOp(point, point, radius^2) -> RBF
        // PointWeightFieldType [index] -> weighting
        template
        <
            class TreeType,
            class RadiusSqrOp,
            class BasisFunctionOp,
            class PointWeightFieldType
        >
        void buildWeightsImpl
        (
            const TreeType& tree,
            const RadiusSqrOp& radiusSqrOp,
            const BasisFunctionOp& basisFuncOp,
            const PointWeightFieldType& pointWeightFld
        );


public:

    //- Runtime type information
    ClassName("filterField");


    // Constructors

        //- Default construct
        FilterField() noexcept = default;

        //- Construct with weights for a field of points
        //- (constant search radius)
        FilterField
        (
            const pointField& points,
            const scalar radius,
            const RBF_type interp = RBF_type::RBF_linear
        );

        //- Construct with weights for faceCentres of a meshedSurface
        //- using a constant search radius search
        //- (or relative radius multiplier for the face bounding spheres).
        FilterField
        (
            const meshedSurface& geom,
            const scalar radius,
            const bool relative = false,
            const RBF_type interp = RBF_type::RBF_linear
        );


    // Member Functions

        //- Reset to unweighted (pass-through)
        void reset();

        //- Create weights for field of points (constant search radius)
        void reset
        (
            const pointField& points,
            const scalar radius,
            const RBF_type interp = RBF_type::RBF_linear
        );

        //- Create weights for meshedSurface using a constant search radius
        //- or optionally with a search radius multiplier for the
        //- face bounding spheres.
        void reset
        (
            const meshedSurface& geom,
            const scalar radius,
            const bool relative = false,
            const RBF_type interp = RBF_type::RBF_linear
        );


    // Evaluation

        //- Return the median smoothed field
        template<class Type>
        tmp<Field<Type>> evaluate
        (
            const tmp<Field<Type>>& tinput,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MappedFileFilterFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
