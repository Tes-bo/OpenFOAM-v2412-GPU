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
    Foam::functionObjects::runTimePostPro::contourFilter

Description
    Iso-surface contours of OpenFOAM volume fields.

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        type        | The surface type: isoSurface          | yes |
        field       | The field defining the surface        | yes |
        colourField | The field to display on the surface   | yes |
        values      | List of iso-values to define the surface(s) | yes |
    \endtable

SourceFiles
    contourFilter.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_contourFilter_H
#define functionObjects_runTimePostPro_contourFilter_H

#include "volumeFilter.H"
#include "fieldVisualisationBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                        Class contourFilter Declaration
\*---------------------------------------------------------------------------*/

class contourFilter
:
    public volumeFilter,
    public fieldVisualisationBase
{
protected:

    // Protected Data

        //- Name of field to colour by
        word colourFieldName_;

        //- The iso values
        List<scalar> values_;


    // Protected Member Functions

        //- No copy construct
        contourFilter(const contourFilter&) = delete;

        //- No copy assignment
        void operator=(const contourFilter&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("isoSurface");


    // Constructors

        //- Construct from dictionary
        contourFilter
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~contourFilter() = default;


    // Member Functions

        //- Add cutting planes to scene (using simulation source)
        bool addGeometry
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Add cutting planes to scene (using simulation source)
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Add cutting planes to scene (using simulation source)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
