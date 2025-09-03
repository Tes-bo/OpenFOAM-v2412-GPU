/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Foam::cellToFaceZone

Description
    A \c topoSetSource to select faces with only one
    neighbour (i.e. outside) in a specified \c cellSet.

    This is just a shortcut for
    - extracting outside faces as a \c faceSet (\c cellToFace with \c outside).
    - using \c setsToFaceZone to convert \c faceSet
    and \c cellSet to oriented \c faceZone.

    Operands:
    \table
      Operand   | Type       | Location
      input     | cellSet(s) | constant/polyMesh/sets/\<set\>
      output    | faceZone   | constant/polyMesh/faceZones
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory entries
        name        <name>;
        type        faceZoneSet;
        action      <action>;
        source      cellToFaceZone;

        // Select one of the below

        // Option-1
        sets
        (
            <word>
            <word>
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
        set         <cellSetName>;

        // Option-4
        zone        <cellZoneName>;

        // Optional entries
        flip        <bool>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Reqd  | Deflt
      name       | Name of faceZone                    | word |  yes  | -
      type       | Type name: faceZoneSet              | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: cellToFaceZone         | word |  yes  | -
      set(s)     | Name of input cellSet(s) containing the slave cells <!--
                                                   --> | word |  yes  | -
      flip       | Flag to select master/slave cells   | bool |  no   | false
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceZone from selected faces
      add      | Add selected faces of a faceZoneSet into this faceZone
      subtract | Remove selected faces of a faceZoneSet from this faceZone
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                    | Type     | Reqd   | Deflt
      sets     | Names of input cellSets        | wordList | choice | -
      zones    | Names of input cellZones       | wordList | cond'l | -
      set      | Name of input cellSet          | word     | choice | -
      zone     | Name of input cellZone         | word     | cond'l | -
    \endverbatim

Notes
  - \c flip=true sets the orientation of faces
    pointing into the \c cellSet, and vice versa.

SourceFiles
    cellToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellToFaceZone_H
#define Foam_cellToFaceZone_H

#include "topoSetFaceZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellSet;

/*---------------------------------------------------------------------------*\
                        Class cellToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class cellToFaceZone
:
    public topoSetFaceZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;

        //- Whether cellSet is slave cells or master cells
        const bool flip_;


    // Private Member Functions

        //- Select outside faces of cellSet
        void selectFaces
        (
            const bitSet& cSet,
            bitSet& selectedFace,
            bitSet& doFlip
        ) const;


public:

    //- Runtime type information
    TypeName("cellToFaceZone");


    // Constructors

        //- Construct from components
        cellToFaceZone
        (
            const polyMesh& mesh,
            const word& cellSetName,
            const bool flip
        );

        //- Construct from dictionary
        cellToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cellToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cellToFaceZone() = default;


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
