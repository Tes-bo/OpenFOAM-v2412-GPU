/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fileFormats::ABAQUSCore

Description
    Core routines used when reading/writing ABAQUS files.

    Face mappings for abaqus deduced from libmesh internals

    Tet4 cells
    \table
        Face       | OpenFOAM | libmesh | abaqus | starcd
        (1 2 3)    | 0        | 2       | 2      | 5
        (0 3 2)    | 1        | 3       | 3      | 4
        (0 1 3)    | 2        | 1       | 1      | 2
        (0 2 1)    | 3        | 0       | 0      | 0
    \endtable

    Pyr5 cells
    \table
        Face       | OpenFOAM | libmesh | abaqus | starcd
        (0 3 2 1)  | 0        | 4       | n/a    | 0
        (0 4 3)    | 1        | 3       | n/a    | 4
        (3 4 2)    | 2        | 2       | n/a    | 3
        (1 2 4)    | 3        | 1       | n/a    | 5
        (0 1 4)    | 4        | 0       | n/a    | 2
    \endtable

    Prism6 cells
    \table
        Face       | OpenFOAM | libmesh | abaqus | starcd
        (0 2 1)    | 0        | 0       | 0      | 0
        (3 4 5)    | 1        | 4       | 1      | 1
        (0 3 5 2)  | 2        | 3       | 4      | 4
        (1 2 5 4)  | 3        | 2       | 3      | 5
        (0 1 4 3)  | 4        | 1       | 2      | 2
    \endtable

    Hex8 cells
    \table
        Face       | OpenFOAM | libmesh | abaqus | starcd
        (0 4 7 3)  | 0        | 4       | 5      | 4
        (1 2 6 5)  | 1        | 2       | 3      | 5
        (0 1 5 4)  | 2        | 1       | 2      | 2
        (3 7 6 2)  | 3        | 3       | 4      | 3
        (0 3 2 1)  | 4        | 0       | 0      | 0
        (4 5 6 7)  | 5        | 5       | 1      | 1
    \endtable

SourceFiles
    ABAQUSCore.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ABAQUSCore_H
#define Foam_ABAQUSCore_H

#include "Fstream.H"
#include "Enum.H"
#include "Map.H"
#include "face.H"
#include "point.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::ABAQUSCore Declaration
\*---------------------------------------------------------------------------*/

class ABAQUSCore
{
public:

    // Public Data, Declarations

        //- Shape-Type - the values are for internal use only!
        enum shapeType : uint8_t
        {
            abaqusUnknownShape = 0,
            abaqusTria  = 0x03,
            abaqusQuad  = 0x04,
            abaqusTet   = 0x84,
            abaqusPyr   = 0x85,
            abaqusPrism = 0x86,
            abaqusHex   = 0x88
        };


    // Public Functions

        //- Classify named element type (eg, S4R) to known/supported
        //- element types.
        //  The input string must be Uppercase!
        static shapeType getElementType(const string& elemTypeName);

        //- The number of points associated with the element type
        inline static int nPoints(shapeType tag)
        {
            return (tag & 0x3F);
        }

        //- True if element type is not unknown/invalid
        inline static bool isValidType(shapeType tag)
        {
            return tag;
        }

        //- True if element type is a 2D shell
        inline static bool isShellType(shapeType tag)
        {
            return (tag & 0x07) && !(tag & 0x80);
        }

        //- True if element type is a 3D element
        inline static bool isSolidType(shapeType tag)
        {
            return (tag & 0x80);
        }

        //- Is a combined (synthetic) face Id?
        inline static bool isEncodedSolidId(const label combinedId)
        {
            return (combinedId < 0);
        }

        //- Combine solid element Id and side Id into synthetic face Id
        inline static label encodeSolidId(const label id, const label side)
        {
            return -(10 * id + side);
        }

        //- Entangle solid element id from synthetic face Id
        inline static label decodeSolidElementId(const label combinedId)
        {
            return
            (
                isEncodedSolidId(combinedId)
              ? (-combinedId / 10)
              : combinedId
            );
        }

        //- Entangle solid side id from synthetic face Id
        //- Synthesize faceId from solid element Id and sideId
        inline static label decodeSolidSideNum(const label combinedId)
        {
            return
            (
                isEncodedSolidId(combinedId)
              ? (-combinedId % 10)
              : 0
            );
        }


protected:

    // Protected Member Functions

        //- From 1-based to 0-based
        inline static void renumber0_elemId(label& elemId)
        {
            if (isEncodedSolidId(elemId))
            {
                // Eg,
                // 1-based (solid 4, side 2) is -42
                // 0-based (solid 3, side 2) is -32
                elemId += 10;
            }
            else
            {
                --elemId;
            }
        }

        //- Face addressing from ABAQUS faces to OpenFOAM faces.
        //  For hex, prism, tet primitive shapes.
        static const Map<labelList>& abaqusToFoamFaceAddr();


    // Protected Classes

        //- Raw reader structure
        struct readHelper
        {
            // General

                //- Additional verbosity
                bool verbose_;


            // Point Handling

                //- Locations of the points (nodes)
                DynamicList<point> points_;

                //- The 1-based abaqus Id for the point (node)
                DynamicList<label> nodeIds_;


            // Element Handling

                //- The element connectivity.
                //  Initially uses the abaqus node Id (1-based)
                //  but remapped to 0-based compact form later.
                DynamicList<labelList> connectivity_;

                //- The 1-based abaqus Id for the element
                DynamicList<label> elemIds_;

                //- The element types
                DynamicList<ABAQUSCore::shapeType> elemTypes_;

                //- The element set ids
                DynamicList<label> elsetIds_;

                //- Mapping of elem set names
                HashTable<label, string> elsetMap_;


            // Constructos

                //- Default construct without verbosity
                explicit readHelper(bool verbosity = false)
                :
                    verbose_(verbosity)
                {}


            // Member Functions

                //- Clear out contents.
                void clear()
                {
                    points_.clear();
                    nodeIds_.clear();

                    connectivity_.clear();
                    elemTypes_.clear();

                    elemIds_.clear();
                    elsetIds_.clear();
                    elsetMap_.clear();
                }

                //- Add a new element set name or return an existing one.
                //  Case-insensitive.
                label addNewElset(const std::string& setName);


                //- Read an abaqus input file
                void read(ISstream& is);

                //- Read entries within a "*Nodes" section.
                //  Appends to points and nodeIds lists.
                //
                //  \return the number of points read
                label readPoints(ISstream& is);

                //- Read entries within an "*Element" section.
                //  If the shape is known/supported, appends to
                //  connectivity, elemType, elemIds lists.
                //
                //  \return the number of elements read
                label readElements
                (
                    ISstream& is,
                    const ABAQUSCore::shapeType shape,
                    const label setId = 0
                );

                //- Read elements within an "*Surface" section.
                //  If the shape is known/supported, appends to
                //  connectivity, elemType, elemIds lists.
                //
                //  \return the number of elements read
                label readSurfaceElements
                (
                    ISstream& is,
                    const label setId = 0
                );


                //- Remove non-shell elements and compact the points
                void purge_solids();

                //- Compact unused points and relabel connectivity
                void compact_nodes();

                //- Renumber elements from 1-based to 0-based
                void renumber_elements_1to0();
        };


    // Constructors

        //- Default construct
        ABAQUSCore() = default;


public:

    // Public Member Functions

        //- Write '*NODE' header and entries to file, optionally with scaling
        //  This is a no-op for an empty list
        static void writePoints
        (
            Ostream& os,
            const UList<point>& points,
            const scalar scaleFactor = 1.0
        );

        //- Calculate face decomposition for non tri/quad faces
        //
        //  \param points the surface points
        //  \param faces  the surface faces
        //  \param decompOffsets begin/end offsets (size+1) into decompFaces
        //  \param decompFaces  List of non-tri/quad decomposed into triangles
        //
        //  \return number of decomposed faces
        static label faceDecomposition
        (
            const UList<point>& points,
            const UList<face>& faces,
            labelList& decompOffsets,
            DynamicList<face>& decompFaces
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
