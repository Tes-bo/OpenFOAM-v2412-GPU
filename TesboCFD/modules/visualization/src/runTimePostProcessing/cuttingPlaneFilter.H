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
    Foam::functionObjects::runTimePostPro::cuttingPlaneFilter

Description
    Cutting planes of OpenFOAM volume fields.

    Example of text object specification:
    \verbatim
    planes
    {
        type            plane;
        planeType       pointAndNormal;

        pointAndNormalDict
        {
                point   (0 0 0);
                normal  (1 0 0);
        }

        offsets         (0 10 20);

        field    T;
    }
    \endverbatim

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        type        | The surface type: plane               | yes |
        planeType   | Selector for plane description        | yes |
        offsets     | Offets of the origin in the normal direction | no | (0)
        field       | The field to display                  | yes |
    \endtable

SourceFiles
    cuttingPlaneFilter.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_cuttingPlaneFilter_H
#define functionObjects_runTimePostPro_cuttingPlaneFilter_H

#include "plane.H"
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
                    Class cuttingPlaneFilter Declaration
\*---------------------------------------------------------------------------*/

class cuttingPlaneFilter
:
    public volumeFilter,
    public fieldVisualisationBase
{
protected:

    // Protected Data

        //- The definition of the plane
        plane plane_;

        //- The offsets to the plane - defaults to (0).
        List<scalar> values_;


    // Protected Member Functions

        //- No copy construct
        cuttingPlaneFilter(const cuttingPlaneFilter&) = delete;

        //- No copy assignment
        void operator=(const cuttingPlaneFilter&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("plane");


    // Constructors

        //- Construct from dictionary
        cuttingPlaneFilter
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~cuttingPlaneFilter() = default;


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

        //- Cleanup files etc.
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
