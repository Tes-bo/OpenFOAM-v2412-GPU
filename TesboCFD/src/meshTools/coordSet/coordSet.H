/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::coordSet

Description
    Holds list of sampling positions

SourceFiles
    coordSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSet_H
#define Foam_coordSet_H

#include "pointField.H"
#include "scalarList.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class coordSet Declaration
\*---------------------------------------------------------------------------*/

class coordSet
:
    public pointField
{
public:

    // Public Data Types

        //- Enumeration defining the output format for coordinates
        enum class coordFormat
        {
            X,          //!< Use 'x' component of points for (scalar) axis
            Y,          //!< Use 'y' component of points for (scalar) axis
            Z,          //!< Use 'z' component of points for (scalar) axis
            RADIUS,     //!< Use mag of points for (scalar) axis
            DISTANCE,   //!< Use additional distance field for (scalar) axis
            XYZ,        //!< Use x,y,z point values
            /// DEFAULT
        };

        //- String representation of coordFormat enum
        static const Enum<coordFormat> coordFormatNames;


protected:

    // Protected Data

        //- Name
        word name_;

        //- Cumulative distance for "distance" write specifier.
        scalarList distance_;

        //- Axis type
        coordFormat axis_;


    // Protected Member Functions

        //- Check for consistent dimensions of points and curve distance
        void checkDimensions() const;


public:

    // Constructors

        /// //- Default construct
        /// coordSet();

        //- Default construct with name and axis type
        coordSet(const word& name, const coordFormat axisType);

        //- Default construct with name and axis type
        coordSet(const word& name, const word& axis);

        //- Copy construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            const List<point>& points,
            const scalarList& dist
        );

        //- Move construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            List<point>&& points,
            scalarList&& dist
        );


    // Member Functions

    // Access

        //- The coord-set name
        const word& name() const noexcept
        {
            return name_;
        }

        //- The sort axis name
        const word& axis() const
        {
            return coordFormatNames[axis_];
        }

        //- Return the points
        const pointField& points() const noexcept
        {
            return static_cast<const pointField&>(*this);
        }

        //- Return the cumulative distance
        const scalarList& distance() const noexcept
        {
            return distance_;
        }

        //- Return the number of points
        label nPoints() const noexcept
        {
            return pointField::size();
        }


    // Edit

        //- Rename the coordinate set
        void rename(const word& newName)
        {
            name_ = newName;
        }

        //- Copy assign new points
        void setPoints(const List<point>& newPoints)
        {
            static_cast<pointField&>(*this) = newPoints;
        }

        //- Move assign new points
        void setPoints(List<point>&& newPoints)
        {
            static_cast<pointField&>(*this) = std::move(newPoints);
        }

        //- Copy assign the cumulative distance
        void setDistance(const scalarList& dist, const bool check=true)
        {
            distance_ = dist;
            if (check) checkDimensions();
        }

        //- Move assign the cumulative distance
        void setDistance(scalarList&& dist, const bool check=true)
        {
            distance_ = std::move(dist);
            if (check) checkDimensions();
        }


    // Output-related

        //- True if axis specification is a vector
        bool hasVectorAxis() const noexcept;

        //- Get coordinate of point according to axis specification.
        //  If axis="distance" is the distance[index]
        scalar scalarCoord(const label index) const;

        //- Get point according to axis="xyz" specification
        const vector& vectorCoord(const label index) const;

        //- Write to stream
        Ostream& write(Ostream& os) const;


    // Other

        //- Gather and sort.
        //  \return (on master) gathered set and overall sort order
        autoPtr<coordSet> gatherSort(labelList& sortOrder) const;


    // Housekeeping

        //- Return the cumulative distance
        const scalarList& curveDist() const noexcept
        {
            return distance_;
        }

        //- Copy assign the cumulative distance
        void setCurveDist(const scalarList& dist)
        {
            setDistance(dist);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
