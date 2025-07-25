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
    Foam::faceToCell

Description
    A \c topoSetCellSource to select all cells based on usage in given
    \c faceSet(s), e.g. select cells that are the \c owner/neighbour/any
    of the faces in a given \c faceSet.

    Operands:
    \table
      Operand   | Type       | Location
      input     | faceSet(s) | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      faceToCell;
        option      <option>;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        sets
        (
            <faceSetName0>
            <faceSetName1>
            ...
        );

        // Option-2
        zones
        (
            <faceZoneName0>
            <faceZoneName1>
            ...
        );

        // Option-3
        set     <faceSetName>;

        // Option-4
        zone    <faceZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: faceToCell             | word |  yes  | -
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
      all       | Cells that are either owner or neighbour of given faces
      any       | Cells that are either owner or neighbour of given faces
      owner     | Cells that are owner of given faces
      neighbour | Cells that are neighbour of given faces
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      sets     | Names of input faceSets    | wordList | cond'l | -
      zones    | Names of input faceZones   | wordList | cond'l | -
      set      | Name of input faceSet      | word     | cond'l | -
      zone     | Name of input faceZone     | word     | cond'l | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource
    - Foam::faceZoneToCell

SourceFiles
    faceToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceToCell_H
#define Foam_faceToCell_H

#include "topoSetCellSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceToCell Declaration
\*---------------------------------------------------------------------------*/

class faceToCell
:
    public topoSetCellSource
{
public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            ANY,
            ALL,
            OWNER,
            NEIGHBOUR
        };


private:

    // Private Data

        static const Enum<faceAction> faceActionNames_;

        //- Add usage string
        static addToUsageTable usage_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;

        //- Option
        faceAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        template<class Selector>
        void combineImpl
        (
            topoSet& set,
            const bool add,
            const Selector& faceLabels
        ) const;

        //- Depending on face to cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("faceToCell");


    // Constructors

        //- Construct from components
        faceToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const faceAction option
        );

        //- Construct from dictionary
        faceToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceToCell() = default;


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
