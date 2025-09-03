/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::PDRobstacle

Description
    Obstacle definitions for PDR

SourceFiles
    PDRobstacle.C
    PDRobstacleIO.C
    PDRobstacleRead.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PDRobstacle_H
#define Foam_PDRobstacle_H

#include "InfoProxy.H"
#include "labelPair.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"
#include "boundBox.H"
#include "DynamicList.H"
#include "pointField.H"
#include "volumeType.H"
#include "memberFunctionSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class boundBox;
class PDRobstacle;

Istream& operator>>(Istream& is, PDRobstacle& obs);
Ostream& operator<<(Ostream& os, const InfoProxy<PDRobstacle>& info);

namespace vtk
{
    class surfaceWriter;
}


/*---------------------------------------------------------------------------*\
                         Class PDRobstacle Declaration
\*---------------------------------------------------------------------------*/

class PDRobstacle
{
public:

    //- Obstacle types (legacy numbering)
    enum legacyTypes
    {
        NONE =  0,          //!< Placeholder
        CUBOID_1  =  1,
        CYLINDER  =  2,
        LOUVER_BLOWOFF = 5,
        LOUVRE_BLOWOFF = 5,
        CUBOID    =  6,
        WALL_BEAM =  7,
        GRATING   =  8,
        OLD_INLET =  9,     //!< ignored (old)
        OLD_BLOWOFF = 10,   //!< ignored (old)
        CIRC_PATCH  = 12,
        RECT_PATCH  = 16,
        DIAG_BEAM   = 22,
        IGNITION    = 41,   //!< ignored (old)
        MESH_PLANE  = 46,
        IGNORE = 200
    };


    // Static Data Members

        //- The max blowoff pressure [bar]
        //  Primarily to catch accidental input in Pa or mbar
        static constexpr int maxBlowoffPressure = 10;


    // Data Members

        //- The group-id
        label groupId;

        //- The obstacle type-id
        int typeId;

        //- The x/y/z orientation (0,1,2)
        direction orient;

        //- Bias for position sorting
        scalar sortBias;

        //- The obstacle location.
        //  Lower corner for boxes, end-centre for cylinders
        point pt;

        //- The obstacle dimensions (for boxes)
        vector span;

        // Accessors for cylinders and diagonal blocks

        inline scalar dia()   const { return span[vector::X]; }
        inline scalar theta() const { return span[vector::Y]; }
        inline scalar len()   const { return span[vector::Z]; }

        inline scalar& dia()   { return span[vector::X]; }
        inline scalar& theta() { return span[vector::Y]; }
        inline scalar& len()   { return span[vector::Z]; }

        union
        {
            scalar wa;
            scalar slat_width;
            scalar blowoff_press;
        };
        union
        {
            scalar wb;
            scalar blowoff_time;
        };
        scalar vbkge;
        scalar xbkge;
        scalar ybkge;
        scalar zbkge;

        union
        {
            int blowoff_type;
            int inlet_dirn;
        };

        string identifier;

public:

    // Constructors

        //- Construct zero-initialized
        PDRobstacle();

        //- Read construct as named dictionary
        explicit PDRobstacle(Istream& is);


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            PDRobstacle,
            read,
            dictionary,
            (
                PDRobstacle& obs,
                const dictionary& dict
            ),
            (obs, dict)
        );


    // Static Member Functions

        //- Read obstacle files and add to the lists
        //  \return the total volume
        static scalar legacyReadFiles
        (
            const fileName& obsFileDir,
            const wordList& obsFileNames,
            const boundBox& meshBb,
            // output
            DynamicList<PDRobstacle>& blocks,
            DynamicList<PDRobstacle>& cylinders
        );

        //- Read obstacle files and set the lists
        //  \return the total volume
        static scalar readFiles
        (
            const fileName& obsFileDir,
            const wordList& obsFileNames,
            const boundBox& meshBb,
            // output
            DynamicList<PDRobstacle>& blocks,
            DynamicList<PDRobstacle>& cylinders
        );


    // Member Functions

        //- Read name / dictionary
        bool read(Istream& is);

        //- Read the 'name' identifier if present
        void readProperties(const dictionary& dict);

        //- Obstacle position accessors
        inline scalar x() const { return pt.x(); }
        inline scalar y() const { return pt.y(); }
        inline scalar z() const { return pt.z(); }
        inline scalar& x() { return pt.x(); }
        inline scalar& y() { return pt.y(); }
        inline scalar& z() { return pt.z(); }


        //- Is obstacle type id cylinder-like?
        inline static bool isCylinder(const label id);

        //- Is obstacle cylinder-like?
        inline bool isCylinder() const;

        //- Reset to a zero obstacle
        void clear();

        //- Scale obstacle dimensions by specified scaling factor
        //  Zero and negative factors are ignored
        void scale(const scalar factor);

        //- Volume of the obstacle
        scalar volume() const;

        //- True if the obstacle is considered to be too small
        bool tooSmall(const scalar minWidth) const;

        //- Set values from single-line, multi-column format.
        //  The only input format, but termed 'legacy' since it may
        //  be replaced in the near future.
        //  \return false if the scanning failed or if the obstacle type
        //      is not supported (or no longer supported)
        bool setFromLegacy
        (
            const int groupTypeId,
            const string& buffer,
            const label lineNo = -1,
            const word& inputFile = word::null
        );

        //- Trim obstacle to ensure it is within the specified bounding box
        //- and return the intersection type.
        //  Returns UNKNOWN for unknown types and invalid bounding boxes
        volumeType trim(const boundBox& bb);

        //- Surface (points, faces) representation
        meshedSurface surface() const;

        //- Add pieces to vtp output
        static label addPieces
        (
            vtk::surfaceWriter& surfWriter,
            const UList<PDRobstacle>& list,
            label pieceId = 0
        );

        //- Generate multi-piece VTK (vtp) file of obstacles
        static void generateVtk
        (
            const fileName& outputDir,
            const UList<PDRobstacle>& obslist,
            const UList<PDRobstacle>& cyllist
        );


    // IOstream Operators

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<PDRobstacle> info() const noexcept
        {
            return *this;
        }

        friend Istream& operator>>(Istream& is, PDRobstacle& obs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators

//- Compare according to x0 position
bool operator<(const PDRobstacle& a, const PDRobstacle& b);

//- For list output, assert that no obstacles are identical
inline bool operator!=(const PDRobstacle& a, const PDRobstacle& b)
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace PDRlegacy
{

/*---------------------------------------------------------------------------*\
                    Class obstacleGrouping Declaration
\*---------------------------------------------------------------------------*/

//- Locations for each instance of an obstacle group.
class obstacleGrouping
:
    public DynamicList<point>
{
    //- Number of obstacles counted
    label nObstacle_;

    //- Number of cylinder-like obstacles counted
    label nCylinder_;


public:

    //- Construct null
    obstacleGrouping()
    :
        nObstacle_(0),
        nCylinder_(0)
    {}

    //- Construct with one location (instance)
    explicit obstacleGrouping(const vector& origin)
    :
        obstacleGrouping()
    {
        append(origin);
    }

    //- Clear obstacle count and locations
    void clear()
    {
        nObstacle_ = 0;
        nCylinder_ = 0;
        DynamicList<point>::clear();
    }

    //- Increment the number of obstacles
    void addObstacle()
    {
        ++nObstacle_;
    }

    //- Increment the number of cylinder-like obstacles
    void addCylinder()
    {
        ++nCylinder_;
    }

    //- The number of obstacles
    label nObstacle() const
    {
        return nObstacle_;
    }

    //- The number of cylinder-like obstacles
    label nCylinder() const
    {
        return nCylinder_;
    }

    //- The number of locations x number of obstacles
    label nTotalObstacle() const
    {
        return size() * nObstacle_;
    }

    //- The number of locations x number of cylinder-like obstacles
    label nTotalCylinder() const
    {
        return size() * nCylinder_;
    }

    //- The number of locations x number of obstacles
    label nTotal() const
    {
        return size() * (nObstacle_ + nCylinder_);
    }

    //- Add a location
    using DynamicList<point>::append;

    //- Add a location
    void append(const scalar x, const scalar y, const scalar z)
    {
        append(point(x, y, z));
    }
};


// Service Functions

//- Read obstacle files, do counting only.
//  \return nObstacle, nCylinder read
labelPair readObstacleFiles
(
    const fileName& obsFileDir,
    const wordList& obsFileNames,
    Map<obstacleGrouping>& groups
);


//- Read obstacle files and add to the lists
//  \return the total volume
scalar readObstacleFiles
(
    const fileName& obsFileDir,
    const wordList& obsFileNames,
    const Map<obstacleGrouping>& groups,
    const boundBox& meshBb,
    // output
    DynamicList<PDRobstacle>& blocks,
    DynamicList<PDRobstacle>& cylinders
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDRlegacy


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators

//- Locations for each instance of an obstacle group.
inline Ostream& operator<<
(
    Ostream& os,
    const PDRlegacy::obstacleGrouping& group
)
{
    os  << token::BEGIN_LIST
        << group.size() << token::SPACE
        << group.nObstacle() << token::SPACE
        << group.nCylinder() << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PDRobstacleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
