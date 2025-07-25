/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::lumpedPointInterpolator

Description
    A simple linear interpolator between two locations, which are
    referenced by index.

    When the interpolator is built for based on searching for nearest and
    next-nearest points, the interpolation will typically cover a (0-0.5) range
    rather than a (0-1) range. The (0.5-1) range implies a flip in the
    relationship of nearest vs. next-nearest.

SourceFiles
    lumpedPointInterpolatorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_lumpedPointInterpolator_H
#define Foam_lumpedPointInterpolator_H

#include "Pair.H"
#include "triFace.H"
#include "scalarList.H"
#include "barycentric2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class lumpedPointInterpolator Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointInterpolator
{
    // Private Data

        //- The id of the nearest point
        label nearest_;

        //- The id of the neighbour point(s)
        label next1_, next2_;

        //- The interpolation weight for the neighbour point(s)
        scalar weight1_, weight2_;


public:

    // Constructors

        //- Default construct, with zero weighting and invalid ids
        inline lumpedPointInterpolator();

        //- Construct with nearest id
        explicit inline lumpedPointInterpolator(const label id);


    //- Destructor
    ~lumpedPointInterpolator() = default;


    // Member Functions

    // Access

        //- True if there is an associated nearest point
        inline bool good() const;

        //- Same as good()
        bool valid() const { return good(); }

        //- The nearest control point, or -1 if !good()
        inline label nearest() const;

        //- The first neighbour control point - identical to next1()
        inline label next() const;

        //- The first neighbour control point, or -1 if !good()
        inline label next1() const;

        //- The second neighbour control point, or -1 if !good()
        inline label next2() const;

        //- The weighting for the nearest point
        inline scalar weight0() const;

        //- The weighting for the first neighbour point,
        //- this also corresponds to the logical location (interval 0-1)
        inline scalar weight1() const;

        //- The weighting for the second neighbour point,
        //- this also corresponds to the logical location (interval 0-1)
        inline scalar weight2() const;


    // Edit

        //- Assign the nearest point, clearing any neighbour
        inline void nearest(const label id);

        //- Assign single neighbour control point and its weight
        inline void next(const label id, const scalar weight);

        //- Assign the neighbour control point and its weight
        inline void next
        (
            const label id,
            const scalar weight,
            const label position    //!< Use 0 or 1 for first/second
        );

        //- Assign all control points and their weights
        //  The triFace points [0,1,2] correspond to [nearest,next1,next2],
        //  respectively
        inline void set(const triFace& ids, const barycentric2D& weights);


    // Evalulate

        //- Linear interpolated value between nearest and next locations
        template<class T>
        inline T interpolate(const UList<T>& input) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lumpedPointInterpolatorI.H"

#endif

// ************************************************************************* //
