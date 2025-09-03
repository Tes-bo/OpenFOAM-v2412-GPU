/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::geometryPatches

Description
    Visualisation of OpenFOAM patches and fields.

    \heading Basic Usage
    \table
        Property    | Description                           | Required | Default
        type        | The surface type: patches             | yes |
        patches     | Patches to display (wordRe list)      | yes |
    \endtable

SourceFiles
    geometryPatches.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_geometryPatches_H
#define functionObjects_runTimePostPro_geometryPatches_H

#include "surface.H"
#include "geometrySurface.H"
#include "fieldVisualisationBase.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                       Class geometryPatches Declaration
\*---------------------------------------------------------------------------*/

class geometryPatches
:
    public geometrySurface,
    public fieldVisualisationBase
{
protected:

    // Protected Data

        //- Requested names of patches to process
        wordRes selectPatches_;

        //- Use cell value on patches instead of patch value itself
        bool nearCellValue_;


    // Protected Member Functions

        //- Gather and convert patches to multi-piece dataset with
        //- vtkPolyData for each patch/processor.
        //  For serial, the pieces are gathered to the master.
        vtkSmartPointer<vtkMultiPieceDataSet>
        gatherPatchPieces(const labelListList& patchIds) const;

        //- Gather and convert patch face centres to multi-piece dataset with
        //- vtkPolyData for each patch/processor.
        //  For serial, the pieces are gathered to the master.
        vtkSmartPointer<vtkMultiPieceDataSet>
        gatherPatchFaceCentres(const labelListList& patchIds) const;


    // Adding Fields - multi-piece

        //- Add patch values.
        //  For nCells == nPoints (eg, only has face centres) add as PointData.
        //  \return 0 on failure to map and nCmpt (eg, 1=scalar, 3=vector)
        //  on success.
        template<class Type>
        int addPatchField
        (
            vtkMultiPieceDataSet* multiPiece,
            const labelListList& patchIds,
            const GeometricField<Type, fvPatchField, volMesh>* fldptr,
            const word& fieldName
        ) const;


        //- No copy construct
        geometryPatches(const geometryPatches&) = delete;

        //- No copy assignment
        void operator=(const geometryPatches&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("patches");


    // Constructors

        //- Construct from dictionary
        geometryPatches
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~geometryPatches() = default;


    // Member Functions

        //- Add geometry surface(s) to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Update actors
        virtual void updateActors(const scalar position);

        //- Clear files used to create the object(s)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "geometryPatchesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
