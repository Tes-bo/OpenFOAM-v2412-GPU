/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::regionToCell

Description
    A \c topoSetCellSource to select cells belonging to a topologically
    connected region (that contains given points). If started inside
    a given \c subCellSet keeps to it; if started outside stays outside.

    Operands:
    \table
      Operand   | Type    | Location
      input     | region  | $FOAM_CASE/constant/\{\<region\>, polyMesh\}
      output    | cellSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        cellSet;
        action      <action>;

        // Mandatory entries
        source      regionToCell;
        insidePoints
        (
            (<p1x> <p1y> <p1z>)
            (<p2x> <p2y> <p2z>)
            ...
        );

        // Optional entries
        set       <cellSetName>;
        zone      <cellZoneName>;
        nErode    <label>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: regionToCell           | word |  yes  | -
      insidePoints  | Coordinate(s) that is inside connected region <!--
                                             --> | vectorList |  yes  | -
      set / zone | Name of cellSet / zone giving mesh subset  | word |  no | -
      nErode     | Number of cell layers to erode mesh to detect holes <!--
           --> in the mesh - set to 0 if not used      | label | no   | 0
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    regionToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionToCell_H
#define Foam_regionToCell_H

#include "topoSetCellSource.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionSplit;

/*---------------------------------------------------------------------------*\
                        Class regionToCell Declaration
\*---------------------------------------------------------------------------*/

class regionToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of cellSet or cellZone to keep to
        word setName_;

        //- Is name a set or a zone
        bool isZone_;

        //- Coordinate(s) that is inside connected region
        const pointField insidePoints_;

        //- Number of layers to erode
        const label nErode_;


    // Private Member Functions

        //- Mark faces inbetween selected and unselected elements
        void markRegionFaces
        (
            const boolList& selectedCell,
            boolList& regionFace
        ) const;

        //- Determine for every disconnected region in the mesh whether
        //- it contains a locationInMesh
        boolList findRegions(const bool verbose, const regionSplit&) const;

        //- Unselect regions not containing a locationInMesh
        void unselectOutsideRegions(boolList& selectedCell) const;

        //- Unselect one layer of cells from selectedCell
        void shrinkRegions(boolList& selectedCell) const;

        //- Erode a given number of layers from selectedCell. Remove any
        //- region that gets disconnected that way.
        void erode(boolList& selectedCell) const;

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("regionToCell");


    // Constructors

        //- Construct from components
        regionToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const pointField& insidePoints,
            const label nErode
        );

        //- Construct from dictionary
        regionToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        regionToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~regionToCell() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
