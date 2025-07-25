/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::boxToCell

Description
    A \c topoSetCellSource to select all cells
    whose cell centre inside given bounding box(es).

    Operands:
    \table
      Operand   | Type    | Location
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
        source      boxToCell;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        boxes
        (
            (<minX> <minY> <minZ>) (<maxX> <maxY> <maxZ>) // box-1
            (<minX> <minY> <minZ>) (<maxX> <maxY> <maxZ>) // box-2
            ...
        );

        // Option-2
        box (<minX> <minY> <minZ>) (<maxX> <maxY> <maxZ>);

        // Option-3
        min (<minX> <minY> <minZ>);
        max (<maxX> <maxY> <maxZ>);

        // Option-4
        min  (<minX> <minY> <minZ>);
        span (<sizeX> <sizeY> <sizeZ>);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: boxToCell              | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                    | Type       | Req'd  | Dflt
      boxes    | Multiple bounding boxes        | vectorList | cond'l | -
      box      | A single bounding box          | vectors    | cond'l | -
      min      | Minimum point (single box)     | vector     | cond'l | -
      max      | Maximum point (single box)     | vector     | cond'l | -
      span     | Box width (single box)         | vector     | cond'l | -
    \endverbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c boxes, \c box, \c min-max or \c min-span pair
    (compatibility with \c searchableBox).

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    boxToCell.C

\*---------------------------------------------------------------------------*/

#ifndef boxToCell_H
#define boxToCell_H

#include "topoSetCellSource.H"
#include "treeBoundBoxList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class boxToCell Declaration
\*---------------------------------------------------------------------------*/

class boxToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Bounding boxes
        treeBoundBoxList bbs_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("boxToCell");


    // Constructors

        //- Construct from components, copying bounding boxes
        boxToCell(const polyMesh& mesh, const treeBoundBoxList& bbs);

        //- Construct from components, moving bounding boxes
        boxToCell(const polyMesh& mesh, treeBoundBoxList&& bbs);

        //- Construct from dictionary
        boxToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from stream - a single box.
        boxToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~boxToCell() = default;


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
