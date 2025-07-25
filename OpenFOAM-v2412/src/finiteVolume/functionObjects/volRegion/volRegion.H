/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::volRegion

Group
    grpFieldFunctionObjects

Description
    Volume (cell) region selection class.

    The adjustments for mesh changes have been implemented with a lazy
    evaluation, to avoid unnecessary recalculation until the values are
    actually required. The update() method is used to ensure the cache
    values are up-to-date.

    Examples of function object specification:
    \verbatim
    volRegion0
    {
        .
        .
        regionType      cellZone;
        name            c0;
        .
        .
    }

    volRegionAll
    {
        .
        .
        regionType      all;
        .
        .
    }
    \endverbatim

Usage
    \table
        Property     | Description                          | Required | Default
        regionType   | Selection type: all/cellSet/cellZone | no | all
        name         | Name of cellSet/cellZone if required | conditional |
    \endtable

See also
    Foam::functionObject

SourceFiles
    volRegion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_volRegion_H
#define Foam_functionObjects_volRegion_H

#include "writeFile.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class polyMesh;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class volRegion Declaration
\*---------------------------------------------------------------------------*/

class volRegion
{
public:

    // Public Data Types

        //- Region type enumeration
        enum regionTypes : char
        {
            vrtAll = 0,         //!< All cells
            vrtCellSet,         //!< A cellSet
            vrtCellZone         //!< A cellZone
        };

        //- Region type names
        static const Enum<regionTypes> regionTypeNames_;


private:

    // Private Member Data

        //- Reference to the volume mesh
        const fvMesh& volMesh_;

        //- The cell ids, from cellSet or multiple cell zones
        labelList cellIds_;

        //- Region IDs (zone ID, ...)
        labelList regionIDs_;

        //- Cached total number of cells selected
        label nCells_;

        //- Cached total selection volume
        scalar V_;

        //- Flag to indicate whether the volRegion requires updating
        bool requireUpdate_;


    // Private Member Functions

        //- Update cellIds, nCells, volume
        void calculateCache();


protected:

    // Protected Data

        //- Region type
        regionTypes regionType_;

        //- Region name (cellSet, cellZone, ...)
        wordRe regionName_;


    // Protected Member Functions

        //- Use all cells, not the cellIDs
        inline bool useAllCells() const noexcept;

        //- Output file header information
        void writeFileHeader(const writeFile& wf, Ostream& file) const;


public:

    //- Run-time type information
    TypeName("volRegion");


    // Constructors

        //- Construct from fvMesh and dictionary
        volRegion(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~volRegion() = default;


    // Member Functions

        //- The region type
        inline regionTypes regionType() const noexcept;

        //- Return the local list of cell IDs.
        //  Empty or invalid for 'all'
        const labelList& cellIDs() const;

        //- Return the total number of cells selected in the region
        inline label nCells() const;

        //- Return total volume of the selected region
        inline scalar V() const;

        //- Update the cached values as required
        //  \return False if the values were already up to date
        bool update();


        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volRegionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
