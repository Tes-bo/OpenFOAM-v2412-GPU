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
    Foam::pointToCell

Description
    A \c topoSetCellSource to select cells with any
    point or any edge within a given \c pointSet(s).

    Operands
    \table
      Operand   | Type         | Location
      input     | pointSeti(s) | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | cellSet      | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      pointToCell;
        option      <option>;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        sets
        (
            <pointSetName1>
            <pointSetName2>
            ...
        );

        // Option-3
        zones
        (
            <pointZoneName0>
            <pointZoneName1>
            ...
        );

        // Option-3
        set     <pointSetName>;

        // Option-4
        zone    <pointZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: pointToCell            | word |  yes  | -
      option     | Selection type - see below          | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Options for the \c option entry:
    \verbatim
      any     | Cells using any point in pointSet
      edge    | Cells using an edge with both points in pointSet
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      sets     | Names of input pointSets   | wordList | cond'l | -
      zones    | Names of input pointZones  | wordList | cond'l | -
      set      | Name of input pointSet     | word     | cond'l | -
      zone     | Name of input pointZone    | word     | cond'l | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    pointToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointToCell_H
#define Foam_pointToCell_H

#include "topoSetCellSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointToCell Declaration
\*---------------------------------------------------------------------------*/

class pointToCell
:
    public topoSetCellSource
{
public:
        //- Enumeration defining the valid options
        enum pointAction
        {
            ANY,    // Cells using any point in set
            EDGE    // Cells using an edge with both points in set
            //ALL   // Possible extension: cells whose all points are in set
        };


private:

    //Private Data

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<pointAction> pointActionNames_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;

        //- Selection type
        pointAction option_;


    // Private Member Functions

        //- Depending on point to cell option add to or delete from cellSet.
        template<class Selector>
        void combineImpl
        (
            topoSet& set,
            const bool add,
            const Selector& pointLabels
        ) const;

        //- Depending on point-to-cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("pointToCell");


    // Constructors

        //- Construct from components
        pointToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const pointAction option
        );

        //- Construct from dictionary
        pointToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        pointToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~pointToCell() = default;


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
