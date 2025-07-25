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
    Foam::functionObjects::runTimePostPro::pointData

Description
    Visualisation of point data

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        representation| sphere/vector                       | yes |
        pointColour   | Override point colour               | no  |
        maxGlyphLength | Limit for glyph representation     | yes |
    \endtable

SourceFiles
    pointData.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_pointData_H
#define functionObjects_runTimePostPro_pointData_H

#include "geometryBase.H"
#include "Enum.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkActor;
class vtkPolyData;
class vtkPolyDataMapper;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                          Class pointData Declaration
\*---------------------------------------------------------------------------*/

class pointData
:
    public geometryBase
{
public:

    // Public Enumerations

        //- Point representation types
        enum representationType
        {
            rtSphere,           //!< "sphere"
            rtVector            //!< "vector"
        };

        //- Names for point representation types
        static const Enum<representationType> representationTypeNames;


protected:

    // Protected Data

        //- Representation type
        representationType representation_;

        //- Max glyph length
        scalar maxGlyphLength_;

        //- Point colour
        autoPtr<Function1<vector>> pointColour_;


    // Protected Member Functions

        //- Add the point data to the renderer
        void addPoints
        (
            const label framei,
            vtkActor* actor,
            vtkPolyDataMapper* mapper,
            vtkPolyData* data
        ) const;

        //- No copy construct
        pointData(const pointData&) = delete;

        //- No copy assignment
        void operator=(const pointData&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("pointData");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointData,
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
        pointData
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    // Selectors

        //- Return selected pointData
        static autoPtr<pointData> New
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours,
            const word& pointDataName
        );


    //- Destructor
    virtual ~pointData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
