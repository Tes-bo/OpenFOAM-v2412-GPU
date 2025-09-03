/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

Namespace
    Foam::functionObjects::runTimePostPro

Description
    Classes and objects used in the implementation of
    Foam::functionObjects::runTimePostProcessing

Class
    Foam::functionObjects::runTimePostProcessing

Group
    grpGraphicsFunctionObjects

Description
    Generate images during run-time.

    The functionality makes use of the VTK libraries (see http://www.vtk.org)
    which provide a broad set of functionality for scene composition and
    manipulation.

    Images are generated using a combination of function object output
    and additional data - eg, surfaces and text.
    Current capabilities include support for:
    - Camera
    - Objects
    - Points
    - Lines
    - Surfaces
    - Scalar bars
    - Annotations
    - Selection of colour maps
    .

    Mandatory entries:
    \verbatim
    output
    {
        name        image;
        width       800;
        height      600;
    }

    camera
    {
        ...
    }

    colours
    {
        ...
    }

    text
    {
    }
    \endverbatim

    Optional entries:
    \verbatim
    points
    {
    }

    lines
    {
        ...
    }

    surfaces
    {
        ...
    }
    \endverbatim

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        type        | Type-name: runTimePostProcessing      | yes |
        debug       | Additional debug information          | no  | false
        parallel    | Allow parallel rendering              | no  | true
    \endtable

SourceFiles
    runTimePostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostProcessing_H
#define functionObjects_runTimePostProcessing_H

#include "fvMeshFunctionObject.H"
#include "objectRegistry.H"
#include "mapPolyMesh.H"
#include "PtrList.H"
#include "scene.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkMultiPieceDataSet;
class vtkMultiProcessController;
class vtkRenderWindow;
class vtkRenderer;
template<class T> class vtkSmartPointer;

namespace Foam
{
namespace functionObjects
{
// Forward Declarations
class runTimePostProcessing;

namespace runTimePostPro
{
// Forward Declarations
class pointData;
class pathline;
class surface;
class text;
}

/*---------------------------------------------------------------------------*\
                   Class runTimePostProcessing Declaration
\*---------------------------------------------------------------------------*/

class runTimePostProcessing
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Output information
        struct outputType
        {
            word name_;
            label width_;
            label height_;
        };

        //- Output instance
        outputType output_;

        //- Allow parallel rendering
        bool parallel_;

        //- Scene manager
        runTimePostPro::scene scene_;

        //- List of points
        PtrList<runTimePostPro::pointData> points_;

        //- List of lines
        PtrList<runTimePostPro::pathline> lines_;

        //- List of surfaces
        PtrList<runTimePostPro::surface> surfaces_;

        //- List of text
        PtrList<runTimePostPro::text> text_;


    // Private Member Functions

        //- Helper function to read scene objects
        template<class Type>
        void readObjects(const dictionary& dict, PtrList<Type>& objects) const;

        //- Construct controller and dispatch to render
        void render();

        //- Static function for SetSingleMethod binding
        static void render(vtkMultiProcessController* ctrl, void* data);

        //- Render scene using given controller
        void render(vtkMultiProcessController* ctrl);


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("runTimePostProcessing");


    // Constructors

        //- Construct from dictionary
        runTimePostProcessing
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~runTimePostProcessing() = default;


    // Member Functions

        //- Allow parallel rendering
        bool parallel() const
        {
            return parallel_;
        }

        //- May need to gather parts to render on single-processor
        //  True when OpenFOAM is running in parallel but VTK is not.
        bool needsCollective() const
        {
            return Pstream::parRun() && !parallel_;
        }

        //- Reference to the underlying OpenFOAM mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }


        //- Read the post-processing controls
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "runTimePostProcessingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
