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

Class
    Foam::functionObjects::runTimePostPro::functionObjectSurface

Description
    Visualisation of surface data from function object output,
    typically the result of a prior sampled surfaces operation.

    Operates in a "live" mode, in which the previously sampled surfaces
    are retrieved from the functionObject registry.
    Or in the traditional file-based mode, in which the sampling is used
    to generate a file and its name is retrieved from the functionObject
    properties.

    File-based import is restricted to "vtk" and "vtp" formats.

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        type        | The type: functionObjectSurface       | yes |
    \endtable

SourceFiles
    functionObjectfunctionObjectSurface.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_functionObjectSurface_H
#define functionObjects_runTimePostPro_functionObjectSurface_H

#include "geometrySurface.H"
#include "functionObjectBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                    Class functionObjectSurface Declaration
\*---------------------------------------------------------------------------*/

class functionObjectSurface
:
    public geometrySurface,
    public functionObjectBase
{
protected:

    // Protected Member Functions

        //- No copy construct
        functionObjectSurface(const functionObjectSurface&) = delete;

        //- No copy assignment
        void operator=(const functionObjectSurface&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("functionObjectSurface");


    // Constructors

        //- Construct from dictionary
        functionObjectSurface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~functionObjectSurface() = default;


    // Member Functions

        //- Add functionObjectSurface to scene (using simulation source)
        bool addGeometry
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Add functionObjectSurface to scene (using file source)
        bool addGeometryFromFile
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Add functionObjectSurface(s) to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Clear files used to create the object(s)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
