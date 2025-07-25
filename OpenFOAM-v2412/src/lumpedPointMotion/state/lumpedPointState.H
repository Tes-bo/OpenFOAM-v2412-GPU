/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::lumpedPointState

Description
    The \a state of lumped points corresponds to positions and rotations.

    This class encapsulates the response from the external application and
    serves as the entry point for applying relaxation, sub-stepping etc.

    \heading Dictionary input format
    \table
        Property    | Description                           | Required | Default
        points      | List of points                        | yes |
        angles      | List of Euler rotation angles         | yes |
        rotationOrder | The Euler-angle rotation order      | no  | zxz
        degrees     | Rotation angles in degrees            | no  | false
    \endtable

    \heading Plain input format.
    Blank and comment lines starting with a '#' character are ignored.
    The angles are always in radians.
    \verbatim
    NumPoints
    x0 y0 z0 eulerz0 eulerx'0 eulerz''0
    x1 y1 z1 eulerz1 eulerx'1 eulerz''1
    ...
    \endverbatim

SeeAlso
    Foam::coordinateRotations::euler, Foam::quaternion

SourceFiles
    lumpedPointState.C
    lumpedPointStateI.H

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointState_H
#define lumpedPointState_H

#include "dictionary.H"
#include "scalarList.H"
#include "pointField.H"
#include "scalarField.H"
#include "vectorField.H"
#include "tensorField.H"
#include "quaternion.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class lumpedPointState Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointState
{
public:

    // Data Types

        //- Input format types
        enum class inputFormatType
        {
            PLAIN,      //!< "plain" is a simple ASCII format
            DICTIONARY  //!< "dictionary" is the OpenFOAM dictionary format
        };


    // Static Data

        //- Names for the input format types
        static const Enum<inputFormatType> formatNames;


private:

    // Private Data

        //- Positions of lumped points
        pointField  points_;

        //- Orientation of lumped points (as Euler angles)
        vectorField angles_;

        //- The Euler-angle rotation order (default: zxz)
        quaternion::eulerOrder order_;

        //- Euler angles measured in degrees
        bool degrees_;

        //- Tensor rotation of lumped points
        mutable std::unique_ptr<tensorField> rotationPtr_;


    // Private Member Functions

        void calcRotations() const;

        void readDict
        (
            const dictionary& dict,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

public:

    // Public Data

        //- Enable/disable visualization of unused points
        static bool visUnused;

        //- The length for visualization triangles
        static scalar visLength;


    // Constructors

        //- Default construct
        lumpedPointState();

        //- Copy construct
        lumpedPointState(const lumpedPointState& rhs);

        //- Copy construct from points and angles
        lumpedPointState
        (
            const pointField& pts,
            const vectorField& ang,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

        //- Copy construct from points with zero-rotation
        explicit lumpedPointState
        (
            const pointField& pts,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

        //- Construct from points with zero-rotation
        explicit lumpedPointState
        (
            tmp<pointField>& pts,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

        //- Construct from dictionary
        explicit lumpedPointState
        (
            const dictionary& dict,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );


    //- Destructor
    virtual ~lumpedPointState() = default;


    // Member Functions

        //- Has positions and consistent number of rotations?
        inline bool good() const;

        //- If no points were specified
        inline bool empty() const;

        //- The number of points
        inline label size() const;

        //- Same as good()
        bool valid() const { return good(); }

        //- The points corresponding to mass centres
        inline const pointField& points() const;

        //- The orientation of the points (mass centres)
        inline const vectorField& angles() const;

        //- The local-to-global transformation for each point
        inline const tensorField& rotations() const;

        //- Scale points by given factor.
        //  Zero and negative values are ignored.
        void scalePoints(const scalar scaleFactor);

        //- The Euler-angle rotation order
        inline quaternion::eulerOrder rotationOrder() const;

        //-  Rotation angles in degrees
        inline bool degrees() const;

        //- Relax the state
        //  alpha = 1 : no relaxation
        //  alpha < 1 : relaxation
        //  alpha = 0 : do nothing
        void relax(const scalar alpha, const lumpedPointState& prev);

        //- Read input as dictionary content
        bool readData
        (
            Istream& is,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

        //- Output as dictionary content
        bool writeData(Ostream& os) const;

        //- Output as dictionary content
        void writeDict(Ostream& os) const;

        //- Read input as plain content
        bool readPlain
        (
            Istream& is,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

        //- Output as plain content
        void writePlain(Ostream& os) const;

        //- Read input from file (master process only) using specified format
        bool readData
        (
            const inputFormatType fmt,
            const fileName& file,
            const quaternion::eulerOrder rotOrder = quaternion::eulerOrder::ZXZ,
            const bool degrees = false
        );

        //- Output points/rotations as VTK file for debugging/visualization
        //  The points are written as vertices, rotation as a triangle
        void writeVTP
        (
            const fileName& outputFile,
            const labelListList& lines = labelListList(),
            const labelList& pointIds = labelList::null()
        ) const;


    // Member Operators

        //- Assignment operator
        void operator=(const lumpedPointState& rhs);

        //- Shift points by specified origin
        void operator+=(const point& origin);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lumpedPointStateI.H"

#endif

// ************************************************************************* //
