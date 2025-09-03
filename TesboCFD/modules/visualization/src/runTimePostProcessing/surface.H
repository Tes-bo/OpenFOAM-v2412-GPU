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
    Foam::functionObjects::runTimePostPro::surface

Description
    Visualisation of surface data with additional routines for handling
    parallel distributed data.

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        representation| none/glyph/wireframe/surface/surfaceWithEdges | yes |
        surfaceColour | Override surface colour             | no  |
        edgeColour    | Override edge colour                | no  |
        featureEdges  | Display surface feature edges       | no  | false
        maxGlyphLength | Limit for glyph representation     | yes | 0
        backFaceCulling | Cull back face                    | no  | false
        frontFaceCulling | Cull front face                  | no  | false
    \endtable

SourceFiles
    surface.C
    surfaceGather.C
    surfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_surface_H
#define functionObjects_runTimePostPro_surface_H

#include "geometryBase.H"
#include "DimensionedField.H"
#include "Enum.H"
#include "runTimeSelectionTables.H"

#include "vtkSmartPointer.h"
#include "vtkMultiPieceDataSet.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkActor;
class vtkRenderer;
class vtkCellData;
class vtkCompositeDataGeometryFilter;
class vtkFeatureEdges;
class vtkPointData;
class vtkPolyData;

namespace Foam
{
// Forward Declarations
class polySurface;
class polySurfaceGeoMesh;
class polySurfacePointGeoMesh;
}


// These need to shift elsewhere

namespace Foam
{
namespace vtk
{
namespace Tools
{

//- Functional call with null-pointer check
vtkCellData* GetCellData(vtkDataSet* dataset);

//- Functional call with null-pointer check
vtkPointData* GetPointData(vtkDataSet* dataset);


//- Default field access is vtkCellData
template<class Type>
struct FieldAccess
{
    vtkCellData* operator()(vtkDataSet* dataset) const
    {
        return Tools::GetCellData(dataset);
    }
};


// Specializations on OpenFOAM type

//- PointAccess for point fields (on polySurfacePointGeoMesh)
template<>
struct FieldAccess<::Foam::polySurfacePointGeoMesh>
{
    vtkPointData* operator()(vtkDataSet* dataset) const
    {
        return Tools::GetPointData(dataset);
    }
};

} // End namespace Tools
} // End namespace vtk
} // End namespace Foam


// More


namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                          Class surface Declaration
\*---------------------------------------------------------------------------*/

class surface
:
    public geometryBase
{
public:

    // Public enumerations

        //- Surface representation types
        enum representationType
        {
            rtNone,             //!< "none"
            rtGlyph,            //!< "glyph"
            rtWireframe,        //!< "wireframe"
            rtSurface,          //!< "surface"
            rtSurfaceWithEdges  //!< "surfaceWithEdges"
        };

        //- Names for surface representation types
        static const Enum<representationType> representationTypeNames;


protected:

    // Protected Data

        //- Representation type
        representationType representation_;

        //- Activate feature edges
        bool featureEdges_;

        //- Back face culling option; default = off
        bool backFaceCulling_;

        //- Front face culling option; default = off
        bool frontFaceCulling_;

        //- Surface colour
        autoPtr<Function1<vector>> surfaceColour_;

        //- Edge colour
        autoPtr<Function1<vector>> edgeColour_;

        //- Surface actor
        vtkSmartPointer<vtkActor> surfaceActor_;

        //- Edge actor
        vtkSmartPointer<vtkActor> edgeActor_;

        //- Max glyph length for representation type rtGlyph
        scalar maxGlyphLength_;


    // Protected Member Functions

        //- Set the representation
        void setRepresentation(vtkActor* actor) const;

        //- Add feature edges to scene
        void addFeatureEdges
        (
            vtkRenderer* renderer,
            vtkFeatureEdges* featureEdges
        ) const;

        //- Add feature edges to scene
        void addFeatureEdges
        (
            vtkRenderer* renderer,
            vtkPolyData* data
        ) const;

        //- Add feature edges to scene
        void addFeatureEdges
        (
            vtkRenderer* renderer,
            vtkCompositeDataGeometryFilter* input
        ) const;


        //- Gather and convert polySurface to multi-piece dataset with
        //- vtkPolyData for the leaves.
        //  If VTK is also running in parallel, each surface is left
        //  as a processor-local piece. Otherwise all processor-local
        //  surfaces are gathered onto the master in their correponding
        //  slots.
        vtkSmartPointer<vtkMultiPieceDataSet>
        gatherSurfacePieces(const polySurface* surf) const;

        //- Gather and convert polySurface to multi-piece dataset with
        //- vtkPolyData for the leaves.
        //  If VTK is also running in parallel, each surface is left
        //  as a processor-local piece. Otherwise all processor-local
        //  surfaces are gathered onto the master in their correponding
        //  slots.
        vtkSmartPointer<vtkMultiPieceDataSet>
        gatherFaceCentres(const polySurface* surf) const;


    // Adding Fields - single-piece

        //- Add field of Type to piece as VTK field data in GeoMeshType slot.
        //  GeoMeshType distinguishes between vtkCellData and vtkPointData
        template<class Type, class GeoMeshType>
        bool addField
        (
            vtkDataSet* piece,       //!< The VTK piece (null protected)
            const Field<Type>& fld,  //!< The field values to add
            const word& fieldName    //!< The field name to use
        ) const;

        //- Attempt cast of regIOobject to DimensionedField\<Type\> and
        //- add to piece as VTK field data in GeoMeshType slot.
        //  GeoMeshType distinguishes between vtkCellData and vtkPointData
        template<class Type, class GeoMeshType>
        bool addDimField
        (
            vtkDataSet* piece,       //!< The VTK piece (null protected)
            const regIOobject* ioptr, //!< The field values to add
            const word& fieldName    //!< The field name to use
        ) const;

        //- Attempt cast of regIOobject to standard DimensionedField types
        //- and add to piece when possible
        template<class GeoMeshType>
        bool addDimField
        (
            vtkDataSet* piece,       //!< The VTK piece (null protected)
            const regIOobject* ioptr, //!< The field values to add
            const word& fieldName    //!< The field name to use
        ) const;


    // Adding Fields - multi-piece

        //- Add DimensionedField of Type to multi-piece as VTK field data in
        //- GeoMeshType slot (CELL | POINT).
        template<class Type, class GeoMeshType>
        bool addDimField
        (
            vtkMultiPieceDataSet* multiPiece,
            const DimensionedField<Type, GeoMeshType>* fldptr,
            const word& fieldName
        ) const;

        //- Attempt cast of regIOobject to DimensionedField\<Type\> and
        //- add in multi-piece as VTK field data in
        //- GeoMeshType slot (CELL | POINT).
        template<class Type, class GeoMeshType>
        bool addDimField
        (
            vtkMultiPieceDataSet* multiPiece, //!< The VTK pieces
            const regIOobject* ioptr, //!< The field values to add
            const word& fieldName    //!< The field name to use
        ) const;

        //- Attempt cast of regIOobject to standard DimensionedField types
        //- and add when possible in GeoMeshType slot (CELL | POINT).
        template<class GeoMeshType>
        bool addDimField
        (
            vtkMultiPieceDataSet* multiPiece, //!< The VTK pieces
            const regIOobject* ioptr, //!< The field values to add
            const word& fieldName    //!< The field name to use
        ) const;

        //- Add using regIOobject information obtained from surface
        template<class GeoMeshType>
        bool addDimField
        (
            vtkMultiPieceDataSet* multiPiece,
            const polySurface* surf,
            const word& fieldName
        ) const;


        //- No copy construct
        surface(const surface&) = delete;

        //- No copy assignment
        void operator=(const surface&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("surface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surface,
            dictionary,
            (
                const runTimePostProcessing& parent,
                const dictionary& dict,
                const HashPtrTable<Function1<vector>>& colours
            ),
            (parent, dict, colours)
        );


    // Constructors

        //- Construct from dictionary
        surface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    // Selectors

        //- Return selected surface
        static autoPtr<surface> New
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours,
            const word& surfaceName
        );


    //- Destructor
    virtual ~surface();


    // Member Functions

        //- Update the actors
        virtual void updateActors(const scalar position);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
