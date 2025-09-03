/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::fieldVisualisationBase

Description
    Base class for scene objects

    Dictionary controls - colour by field
    \table
        Property    | Description                           | Required | Default
        colourBy    | Colouring type (color / field)        | yes |
        range       | Lower/upper range to display          | yes |
        smooth      | Request smoother output               | no  | false
        colourMap   | Colour map for rendering              | no  | rainbow
        nColours    | Number of colours in colour map       | no  | 256
        scalarBar   | Scalar-bar sub-dictionary             | yes |
    \endtable

Colour maps include "coolToWarm" ("blueWhiteRed"), "coldAndHot",
"fire", "rainbow", "greyscale" ("grayscale"), "xray". For historical
reasons, the default is still "rainbow".

SourceFiles
    fieldVisualisationBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldVisualisationBase_H
#define functionObjects_fieldVisualisationBase_H

#include "dictionary.H"
#include "Tuple2.H"
#include "Enum.H"
#include "vector.H"
#include "MinMax.H"
#include "HashPtrTable.H"
#include "scalarBar.H"
#include "Function1.H"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkActor;
class vtkCompositeDataSet;
class vtkCompositeDataGeometryFilter;
class vtkCompositePolyDataMapper;
class vtkDataSet;
class vtkDataSetAlgorithm;
class vtkFieldData;
class vtkLookupTable;
class vtkMapper;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkRenderer;


namespace Foam
{
namespace functionObjects
{
// Forward Declarations
class runTimePostProcessing;

namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                   Class fieldVisualisationBase Declaration
\*---------------------------------------------------------------------------*/

class fieldVisualisationBase
{
public:

    // Public Enumerations

        //- Colouring type
        enum colourByType : uint8_t
        {
            cbColour,           //!< "colour" : Use specified colour
            cbField             //!< "field" : Use named field
        };

        //- Enumeration names for colourByType
        static const Enum<colourByType> colourByTypeNames;

        //- Colour map enumerations
        enum colourMapType : uint8_t
        {
            cmCoolToWarm,       //!< ParaView "Cool To Warm" blue-white-read
            cmBlueWhiteRed = cmCoolToWarm,
            cmColdAndHot,       //!< ParaView "Cold and Hot"
            cmFire,             //!< ParaView "Black-Body Radiation"
            cmRainbow,          //!< "rainbow"
            cmGreyscale,        //!< ParaView "Grayscale"
            cmXray              //!< ParaView "X Ray"
        };

        //- Enumeration names for colourMapType
        static const Enum<colourMapType> colourMapTypeNames;


        //- Enumeration of the data field associations
        //  These values are used internally and do NOT correspond to the
        //  vtkDataObject::FieldAssociations enumeration.
        enum FieldAssociation
        {
            NO_DATA  = 0,           //!< No associated data
            CELL_DATA = 0x1,        //!< Associated with cells (faces)
            POINT_DATA = 0x2,       //!< Associated with points
            CELL_POINT_DATA = 0x3   //!< Associated with cells and/or points
        };


        //- General field characteristics.
        //  For convenience, the interface is exposed but external use is
        //  highly discouraged.
        struct fieldSummary
        {
            int nComponents_;
            unsigned association_;
            scalarMinMax range_;

            //- Default construct
            fieldSummary()
            :
                nComponents_(0),
                association_(0u),
                range_()
            {}

            //- Parallel reduction. A no-op if Pstream::parRun() is false
            void reduce();

            //- True if nComponents_ == 1
            bool isScalar() const
            {
                return nComponents_ == 1;
            }

            //- True if nComponents_ == 3
            bool isVector() const
            {
                return nComponents_ == 3;
            }

            //- True if association_ is non-zero
            bool exists() const
            {
                return association_;
            }

            //- True if there is a POINT_DATA association
            bool hasPointData() const
            {
                return (association_ & FieldAssociation::POINT_DATA);
            }

            InfoProxy<fieldSummary> info() const
            {
                return InfoProxy<fieldSummary>(*this);
            }
        };


protected:

    // Protected Data

        //- Colours
        const HashPtrTable<Function1<vector>>& colours_;

        //- Field name
        word fieldName_;

        //- Requested smoother fields (eg, interpolate cell -> point values)
        bool smooth_;

        //- Colour by type
        colourByType colourBy_;

        //- Colour map type
        colourMapType colourMap_;

        //- Number of colours in colour map
        int nColours_;

        //- Range of values
        Tuple2<scalar, scalar> range_;

        //- Scalar bar characteristics
        scalarBar scalarBar_;


    // Protected Member Functions

        //- Query DataSet for field name and its field association
        static fieldSummary queryFieldSummary
        (
            const word& fieldName,
            vtkDataSet* dataset
        );

        //- Query composite DataSet for field name and its FieldAssociation
        static fieldSummary queryFieldSummary
        (
            const word& fieldName,
            vtkCompositeDataSet* data
        );

        //- Query DataSet for field name and its field association
        static FieldAssociation queryFieldAssociation
        (
            const word& fieldName,
            vtkDataSet* dataset
        );

        //- Query composite DataSet for field name and its FieldAssociation
        static FieldAssociation queryFieldAssociation
        (
            const word& fieldName,
            vtkCompositeDataSet* data
        );


        //- Add "mag(..)" field for filters that only accept scalars
        static void addMagField
        (
            const word& fieldName,
            vtkFieldData* fieldData
        );

        //- Add "mag(..)" field for filters that only accept scalars
        static void addMagField
        (
            const word& fieldName,
            vtkDataSet* dataset
        );

        //- Add "mag(..)" field for filters that only accept scalars
        static void addMagField
        (
            const word& fieldName,
            vtkCompositeDataSet* data
        );


        //- Set the colour map
        void setColourMap(vtkLookupTable* lut) const;

        //- Add scalar bar (if visible) to renderer
        void addScalarBar
        (
            const scalar position,
            vtkRenderer* renderer,
            vtkLookupTable* lut
        ) const;

        //- Set field/configure mapper, add scalar bar
        void setField
        (
            const scalar position,
            const word& colourFieldName,
            const FieldAssociation fieldAssociation,
            vtkMapper* mapper,
            vtkRenderer* renderer
        ) const;

        //- Add glyphs
        void addGlyphs
        (
            const scalar position,
            const word& scaleFieldName,
            const fieldSummary& scaleFieldInfo,
            const word& colourFieldName,
            const fieldSummary& colourFieldInfo,
            const scalar maxGlyphLength,

            vtkPolyData* data,
            vtkActor* actor,
            vtkRenderer* renderer
        ) const;


        //- No copy construct
        fieldVisualisationBase(const fieldVisualisationBase&) = delete;

        //- No copy assignment
        void operator=(const fieldVisualisationBase&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        fieldVisualisationBase
        (
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~fieldVisualisationBase();


    // Member Functions

    // Access

        //- Return the colours
        const HashPtrTable<Function1<vector>>& colours() const;

        //- Return the field name
        const word& fieldName() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects


// Ostream
Ostream& operator<<
(
    Ostream& os,
    const InfoProxy
    <
        functionObjects::runTimePostPro::fieldVisualisationBase::fieldSummary
    >& proxy
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
