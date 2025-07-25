/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::cellToCell

Description
    A \c topoSetCellSource to select all the cells from given \c cellSet(s)
    or \c cellZone(s).

    Operands:
    \table
      Operand   | Type       | Location
      input     | cellSet(s) | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | cellSet    | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      cellToCell;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        sets
        (
            <cellSetName1>
            <cellSetName2>
            ...
        );

        // Option-2
        zones
        (
            <cellZoneName0>
            <cellZoneName1>
            ...
        );

        // Option-3
        set     <cellSetName>;

        // Option-4
        zone    <cellZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: cellToCell             | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                    | Type     | Req'd  | Dflt
      sets     | Names of input cellSets        | wordList | cond'l | -
      zones    | Names of input cellZones       | wordList | cond'l | -
      set      | Name of input cellSet          | word     | cond'l | -
      zone     | Name of input cellZone         | word     | cond'l | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    cellToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellToCell_H
#define Foam_cellToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellToCell Declaration
\*---------------------------------------------------------------------------*/

class cellToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;


public:

    //- Runtime type information
    TypeName("cellToCell");


    // Constructors

        //- Construct from components
        cellToCell(const polyMesh& mesh, const word& setName);

        //- Construct from dictionary
        cellToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cellToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cellToCell() = default;


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
