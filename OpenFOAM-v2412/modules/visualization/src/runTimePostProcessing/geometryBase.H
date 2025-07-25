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
    Foam::functionObjects::runTimePostPro::geometryBase

Description
    Base class for surface, text handling etc.

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        visible     | Display the object                    | no  | yes
        renderMode  | Shading (flat/gouraud/phong)          | no  | gouraud
        opacity     | Object opacity                        | no  | 1.0
        parallel    | Allow parallel rendering              | no  | true
    \endtable

SourceFiles
    geometryBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_geometryBase_H
#define functionObjects_runTimePostPro_geometryBase_H

#include "dictionary.H"
#include "vector.H"
#include "Function1.H"
#include "HashPtrTable.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations (VTK)
class vtkRenderer;
class vtkActor;

namespace Foam
{
namespace functionObjects
{

// Forward Declarations
class runTimePostProcessing;

namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                        Class geometryBase Declaration
\*---------------------------------------------------------------------------*/

class geometryBase
{
public:

    // Public Enumerations

        //- Surface shading types
        enum renderModeType
        {
            rmFlat,             //!< Flat shading
            rmGouraud,          //!< Gouraud shading
            rmPhong             //!< Phong shading
        };

        //- Names for surface shading types
        static const Enum<renderModeType> renderModeTypeNames;


protected:

    // Protected Data

        //- Reference to the parent function object
        const runTimePostProcessing& parent_;

        //- The surface name
        word name_;

        //- Visible flag
        bool visible_;

        //- Allow parallel rendering
        bool parallel_;

        //- Render mode
        renderModeType renderMode_;

        //- Opacity
        autoPtr<Function1<scalar>> opacity_;

        //- Reference to the colours
        const HashPtrTable<Function1<vector>>& colours_;


    // Protected Functions

        //- Initialise actor
        void initialiseActor(vtkActor* actor) const;

        //- No copy construct
        geometryBase(const geometryBase&) = delete;

        //- No copy assignment
        void operator=(const geometryBase&) = delete;


public:

    //- Debug switch
    static int debug;

    // Constructors

        //- Construct from dictionary
        geometryBase
        (
            const runTimePostProcessing& parent_,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~geometryBase();


    // Member Functions

        //- Return the reference to the parent function object
        const runTimePostProcessing& parent() const;

        //- May need to gather geometry parts to render on single-processor
        //  True when OpenFOAM is running in parallel but VTK is not.
        bool needsCollective() const;

        //- Return the name
        const word& name() const;

        //- Return the visible flag
        bool visible() const
        {
            return visible_;
        }

        //- Allow parallel rendering
        bool parallel() const
        {
            return parallel_;
        }

        //- Return the opacity
        scalar opacity(const scalar position) const;

        //- Return reference to the colours
        const HashPtrTable<Function1<vector>>& colours() const;


    // Scene Interaction

        //- Add geometry to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        ) = 0;

        //- Update the actors
        virtual void updateActors(const scalar position) = 0;

        //- Clear any files used to create the object(s)
        virtual bool clear() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
