/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd
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
    Foam::functionObjects::runTimePostPro::pathline

Description
    Visualisation of line data (eg, streamlines)

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        representation| none/line/tube/vector               | yes |
        lineColour    | Override line colour                | no  |
        tubeRadius    | Radius for tube representation      | yes |
    \endtable

SourceFiles
    pathline.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_pathline_H
#define functionObjects_runTimePostPro_pathline_H

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
                          Class pathline Declaration
\*---------------------------------------------------------------------------*/

class pathline
:
    public geometryBase
{
public:

    // Public Enumerations

        //- Line representations
        enum representationType
        {
            rtNone,     //!< "none"
            rtLine,     //!< "line"
            rtTube,     //!< "tube"
            rtVector    //!< "vector"
        };

        //- Names for line representations
        static const Enum<representationType> representationTypeNames;


protected:

    // Protected Data

        //- Representation type
        representationType representation_;

        //- Radius for the tube filter
        scalar tubeRadius_;

        //- Line colour
        autoPtr<Function1<vector>> lineColour_;


    // Protected Member Functions

        //- Add the pathlines to the renderer
        void addLines
        (
            const label framei,
            vtkActor* actor,
            vtkPolyData* data
        ) const;

        //- No copy construct
        pathline(const pathline&) = delete;

        //- No copy assignment
        void operator=(const pathline&) = delete;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("pathline");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pathline,
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
        pathline
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    // Selectors

        //- Return selected pathline
        static autoPtr<pathline> New
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours,
            const word& pathlineName
        );


    //- Destructor
    virtual ~pathline();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
