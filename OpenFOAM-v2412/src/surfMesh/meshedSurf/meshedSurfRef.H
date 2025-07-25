/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::meshedSurfRef

Description
    Implements a meshed surface by referencing another meshed surface
    or faces/points components.

    In addition to the referencing, supports simple moving/scaling
    of points (uses a deep-copy).

\*---------------------------------------------------------------------------*/

#ifndef Foam_meshedSurfRef_H
#define Foam_meshedSurfRef_H

#include "meshedSurf.H"
#include "refPtr.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshedSurfRef Declaration
\*---------------------------------------------------------------------------*/

class meshedSurfRef
:
    public meshedSurf
{
    // Private Data

        //- An external surface reference
        refPtr<meshedSurf> surf_;

        //- Components
        std::reference_wrapper<const pointField> points_;
        std::reference_wrapper<const faceList> faces_;
        std::reference_wrapper<const labelList> zoneIds_;
        std::reference_wrapper<const labelList> faceIds_;

        //- Locations of moved/scaled points (if any)
        pointField newPoints_;


public:

    // Constructors

        //- Default construct
        meshedSurfRef()
        :
            points_(std::cref<pointField>(pointField::null())),
            faces_(std::cref<faceList>(faceList::null())),
            zoneIds_(std::cref<labelList>(labelList::null())),
            faceIds_(std::cref<labelList>(labelList::null()))
        {}

        //- Construct as reference to a meshedSurf
        explicit meshedSurfRef(const meshedSurf& s)
        :
            meshedSurfRef()
        {
            surf_.cref(s);
        }

        //- Construct from components
        meshedSurfRef
        (
            const pointField& points,
            const faceList& faces,
            const labelList& zoneIds = labelList::null(),
            const labelList& faceIds = labelList::null()
        )
        :
            points_(std::cref<pointField>(points)),
            faces_(std::cref<faceList>(faces)),
            zoneIds_(std::cref<labelList>(zoneIds)),
            faceIds_(std::cref<labelList>(faceIds))
        {}


    //- Destructor
    virtual ~meshedSurfRef() = default;


    // Member Functions

        //- Contains a valid reference?
        bool good() const
        {
            return (surf_ || notNull(points_.get()));
        }

        //- Contains a valid reference?
        bool valid() const { return good(); }

        //- The original points used for the surface
        const pointField& points0() const
        {
            return (surf_ ? surf_.cref().points() : points_.get());
        }

        //- The points used for the surface
        virtual const pointField& points() const
        {
            return (newPoints_.empty() ? points0() : newPoints_);
        }

        //- The faces used for the surface
        virtual const faceList& faces() const
        {
            return (surf_ ? surf_.cref().faces() : faces_.get());
        }

        //- Per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return (surf_ ? surf_.cref().zoneIds() : zoneIds_.get());
        }

        //- Per-face identifier (eg, element Id)
        virtual const labelList& faceIds() const
        {
            return (surf_ ? surf_.cref().faceIds() : faceIds_.get());
        }

        //- Invalid by redirecting to null objects
        void clear()
        {
            surf_.reset();
            points_ = std::cref<pointField>(pointField::null());
            faces_ = std::cref<faceList>(faceList::null());
            zoneIds_ = std::cref<labelList>(labelList::null());
            faceIds_ = std::cref<labelList>(labelList::null());
            newPoints_.clear();
        }

        //- Reset surface
        void reset(const meshedSurf& s)
        {
            clear();
            surf_.cref(s);
        }

        //- Reset components
        void reset
        (
            const pointField& points,
            const faceList& faces,
            const labelList& zoneIds = labelList::null(),
            const labelList& faceIds = labelList::null()
        )
        {
            surf_.reset();
            points_ = std::cref<pointField>(points);
            faces_ = std::cref<faceList>(faces);
            zoneIds_ = std::cref<labelList>(zoneIds);
            faceIds_ = std::cref<labelList>(faceIds);
            newPoints_.clear();
        }


        //- Reset changes in point positions
        void resetPoints()
        {
            newPoints_.clear();
        }

        //- Change point positions
        void movePoints(pointField&& pts)
        {
            newPoints_.transfer(pts);
        }

        //- Change point positions
        void movePoints(const tmp<pointField>& tpts)
        {
            newPoints_.clear();
            if (tpts)
            {
                newPoints_ = tpts;
            }
            tpts.clear();
        }

        //- Scale points: ignore unity and non-positive factors
        void scalePoints(const scalar scaleFactor)
        {
            if (scaleFactor > SMALL && !equal(scaleFactor, 1))
            {
                if (newPoints_.empty())
                {
                    newPoints_ = points0();
                }
                newPoints_ *= scaleFactor;
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
