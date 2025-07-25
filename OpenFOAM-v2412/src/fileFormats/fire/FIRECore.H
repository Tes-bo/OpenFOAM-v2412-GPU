/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fileFormats::FIRECore

Description
    Core routines used when reading/writing AVL/FIRE files.

SourceFiles
    FIRECore.C

\*---------------------------------------------------------------------------*/

#ifndef FIRECore_H
#define FIRECore_H

#include "point.H"
#include "labelList.H"
#include "pointField.H"
#include "IOstreams.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// forward declarations
class polyMesh;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class fileFormats::FIRECore Declaration
\*---------------------------------------------------------------------------*/

class FIRECore
{
public:

    // Public Data, Declarations

        //- Selection Types
        enum selectionType
        {
            cellSelection = 2,
            faceSelection = 3
        };

        //- Shape-Type for FIRE (FLMA) files
        enum shapeType
        {
            fireLine  = 1,
            fireTri   = 2,
            fireQuad  = 3,
            fireTet   = 4,
            fireHex   = 5,
            firePyr   = 6,
            firePrism = 8
        };

        //- Enumeration defining the file extensions for 3D types
        enum fileExt3d
        {
            POLY_ASCII,
            POLY_BINARY,
            POLY_ASCII_Z,
            POLY_BINARY_Z
        };


        //- Integer type (binary format)
        typedef int32_t fireInt_t;

        //- Float type (binary format)
        typedef double  fireReal_t;


protected:

    // Protected Data

        static const Enum<fileExt3d>  file3dExtensions;


    // Protected Member Functions

        //- Construct null
        FIRECore() = default;


        //- Read points.
        // This is the first thing to do when reading FPMA,FPMB,FLMA files.
        // Return the number of points read.
        //
        // The file format is as follows:
        // \verbatim
        // NUMBER_OF_VERTICES
        // x0 y0 z0 x1 y1 z1 ... xN-1 yN-1 zN-1
        // \endverbatim
        static label readPoints(ISstream&, pointField&);

public:

    // Public Member Functions

        //- Resolve base file-name for the given file-type
        static fileName fireFileName
        (
            const fileName& baseName,
            const enum fileExt3d
        );


        //- Get an integer (ascii or binary)
        static label getFireLabel(ISstream&);

        //- Get an point x/y/z (ascii or binary)
        static point getFirePoint(ISstream&);

        //- Extract a string (ascii or binary)
        static std::string getFireString(ISstream&);


        //- Write an integer (ascii or binary)
        static void putFireLabel(OSstream&, const label);

        //- Write multiple integers (ascii or binary)
        static void putFireLabels(OSstream&, const labelUList&);

        //- Write an on-the-fly list of integers (ascii or binary)
        static void putFireLabels
        (
            OSstream&,
            const label count,
            const label start
        );


        //- Write a point x/y/z (ascii or binary)
        static void putFirePoint(OSstream&, const point&);

        //- Write a string (ascii or binary)
        static void putFireString(OSstream&, const std::string&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
