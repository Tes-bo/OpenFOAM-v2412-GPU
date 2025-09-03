/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::scene

Description
    Class to control scene construction and provide main rendering loop

Usage
    \verbatim

    // Optional axes dictionary
    axes
    {
        visible         yes;
        colour          (0 0 0);
        fontColour      (0.7 0.5 0.3);
        fontSize        30;
        lengthFactor    0.2;
    }

    camera
    {
        // Total number of frames to generate
        nFrameTotal 1;

        // Parallel projection flag
        parallelProjection no;

        focalPoint  (0 0 0);
        up          (0 1 0);
        position    (0 0 1);

        // Optional entries
        clipBox     (-0.0206 -0.0254 -0.0005) (0.29 0.0254 0.0005);
        viewAngle   20;
        zoom        1.1;
    }

    colours
    {
    }
    \endverbatim

SourceFiles
    scene.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_scene_H
#define functionObjects_runTimePostPro_scene_H

// OpenFOAM includes
#include "dictionary.H"
#include "Function1.H"
#include "vector.H"
#include "point.H"
#include "boundBox.H"
#include "Enum.H"
#include "HashPtrTable.H"

// VTK includes
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkActor;
class vtkAxesActor;
class vtkRenderWindow;
class vtkRenderer;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                           Class scene Declaration
\*---------------------------------------------------------------------------*/

class scene
{
    // Private Data

        //- Reference to the object registry
        const objectRegistry& obr_;

        //- Object name
        const word name_;


    // Private Member Functions

        //- Read camera properties
        void readCamera(const dictionary& dict);

        //- Read colour properties
        void readColours(const dictionary& dict);

        //- Read axes properties
        void readAxes(const dictionary& dict);

        //- Set visibility of all actors on/off
        void setActorVisibility
        (
            vtkRenderer* renderer,
            const bool visible
        ) const;

        //- No copy construct
        scene(const scene&) = delete;

        //- No copy assignment
        void operator=(const scene&) = delete;


protected:

    // Helper struct to manage axes properties
    struct axesProperties
    {
        bool visible;
        vector colour;
        vector fontColour;
        scalar fontSize;
        scalar lengthFactor;

        //- Axes actor
        vtkSmartPointer<vtkAxesActor> actor;
    };


    // Protected Data

        //- Colours
        HashPtrTable<Function1<vector>> colours_;

        // Camera settings

            //- Position
            autoPtr<Function1<point>> cameraPosition_;

            //- Focal point
            autoPtr<Function1<point>> cameraFocalPoint_;

            //- Up direction
            autoPtr<Function1<vector>> cameraUp_;

            //- View angle
            autoPtr<Function1<scalar>> cameraViewAngle_;

            //- Zoom: 1 = do nothing, >1 = zoom in, <1 = zoom out
            //  - perspective mode: reduces view angle
            //  - parallel mode: manipulate parallel scale
            autoPtr<Function1<scalar>> cameraZoom_;


        // Scene management

            //- Clipping box
            boundBox clipBox_;

            //- Clipping box actor
            vtkSmartPointer<vtkActor> clipBoxActor_;

            //- Axes properties;
            axesProperties axes_;

            //- Parallel projection flag
            bool parallelProjection_;

            //- Number of frames
            label nFrameTotal_;

            //- Start position [0-1]
            scalar startPosition_;

            //- Position [0-1]
            scalar position_;

            //- Change in position per frame
            scalar dPosition_;

            //- Index of current frame
            label currentFrameI_;

            //- Name prefix of output
            word outputName_;


    // Protected Member Functions

        void setCamera(vtkRenderer* renderer) const;

        string frameIndexStr() const;


public:

    // Constructors

        //- Construct from components
        scene(const objectRegistry& obr, const word& name);


    //- Destructor
    virtual ~scene();


    // Member Functions

    // Access

        //- Return the colours
        const HashPtrTable<Function1<vector>>& colours() const;

        //- Return the current frame index
        label frameIndex() const;

        //- Return the current position (in range 0-1)
        scalar position() const;


        void read(const dictionary& dict);

        void initialise(vtkRenderer* renderer, const word& outputName);

        //- Main control loop
        bool loop(vtkRenderer* renderer);

        //- Save image to file
        void saveImage(vtkRenderWindow* renderWindow) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
