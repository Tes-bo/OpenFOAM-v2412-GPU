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
    Foam::cellToFace

Description
    A \c topoSetFaceSource to select all the faces from given \c cellSet(s).

    Operands:
    \table
      Operand   | Type       | Location
      input     | cellSet(s) | constant/polyMesh/sets/\<set\>
      output    | faceSet    | constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceSet;
        action      <action>;

        // Mandatory entries
        source      cellToFace;
        option      <option>;

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
      Property   | Description                         | Type |  Reqd | Deflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: cellToFace             | word |  yes  | -
      option     | Selection type - see below          | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected cells of cellSet(s)
      add      | Add selected faces of cellSet(s) into this faceSet
      subtract | Remove selected faces of cellSet(s) from this faceSet
    \endverbatim

    Options for the \c option entry:
    \verbatim
      all     | All faces of cells in the cellSet
      both    | Faces where both neighbours are in the cellSet
      outside | Faces with only one neighbour in the cellSet
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                    | Type     | Reqd   | Deflt
      sets     | Names of input cellSets        | wordList | choice | -
      zones    | Names of input cellZones       | wordList | cond'l | -
      set      | Name of input cellSet          | word     | choice | -
      zone     | Name of input cellZone         | word     | cond'l | -
    \endverbatim

Note
  - The \c outside option applies to the cellSets individually.

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    cellToFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellToFace_H
#define Foam_cellToFace_H

#include "topoSetFaceSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellToFace Declaration
\*---------------------------------------------------------------------------*/

class cellToFace
:
    public topoSetFaceSource
{
public:
        //- Enumeration defining the valid options
        enum cellAction
        {
            ALL,
            BOTH,
            OUTSIDE
        };


private:

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<cellAction> cellActionNames_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;

        //- Selection type
        cellAction option_;


    // Private Member Functions

        //- Depending on face to cell option add to or delete from cellSet.
        template<class Selector>
        void combineImpl
        (
            topoSet& set,
            const bool add,
            const Selector& cellLabels
        ) const;

        //- Depending on face to cell option add to or delete from set.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("cellToFace");


    // Constructors

        //- Construct from components
        cellToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const cellAction option
        );

        //- Construct from dictionary
        cellToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cellToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cellToFace() = default;


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
