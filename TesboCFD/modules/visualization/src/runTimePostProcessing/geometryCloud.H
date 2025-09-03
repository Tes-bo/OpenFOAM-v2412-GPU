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
    Foam::functionObjects::runTimePostPro::geometryCloud

Description
    Visualisation of cloud data from function object output (file-based only).

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        type        | The point type: geometryCloud         | yes |
        cloud         | The cloud name                      | no  |
        field         | The field for glyphs scaling        | no  |
        colourField   | The field to display                | no  |
    \endtable

SourceFiles
    geometryCloud.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_geometryCloud_H
#define functionObjects_runTimePostPro_geometryCloud_H

#include "pointData.H"
#include "functionObjectBase.H"
#include "IOField.H"

#include "vtkSmartPointer.h"
#include "vtkMultiPieceDataSet.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                     Class geometryCloud Declaration
\*---------------------------------------------------------------------------*/

class geometryCloud
:
    public pointData,
    public fieldVisualisationBase
{
protected:

    // Protected Data

        //- Name of geometryCloud
        word cloudName_;

        //- Name of field to colour by
        word colourFieldName_;

        //- Actor
        vtkSmartPointer<vtkActor> actor_;


    // Protected Member Functions

        //- Gather and convert cloud positions with vtkPolyData for the leaves.
        //  If VTK is also running in parallel, each cloud is left
        //  as a processor-local piece. Otherwise all processor-local
        //  parts are gathered onto the master in their correponding
        //  slots.
        vtkSmartPointer<vtkMultiPieceDataSet>
        gatherCloud(const objectRegistry& obrTmp) const;


        //- Add field
        template<class Type>
        bool addField
        (
            vtkDataSet* piece,
            const Field<Type>& fld,
            const word& fieldName
        ) const;

        //- Add field
        template<class Type>
        bool addCloudField
        (
            vtkMultiPieceDataSet* multiPiece,
            const IOField<Type>* fldptr,
            const word& fieldName
        ) const;

        //- Add field
        template<class Type>
        bool addCloudField
        (
            vtkMultiPieceDataSet* multiPiece,
            const regIOobject* fieldPtr,
            const word& fieldName
        ) const;

        //- Add field
        bool addCloudField
        (
            vtkMultiPieceDataSet* multiPiece,
            const objectRegistry& obrTmp,
            const word& fieldName
        ) const;


        //- No copy construct
        geometryCloud(const geometryCloud&) = delete;

        //- No copy assignment
        void operator=(const geometryCloud&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("geometryCloud");


    // Constructors

        //- Construct from dictionary
        geometryCloud
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~geometryCloud();


    // Member Functions

        //- Add cloud to scene (from simulation)
        bool addGeometry
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Add cloud to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Update actors
        virtual void updateActors(const scalar position);

        //- No-op
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "geometryCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
