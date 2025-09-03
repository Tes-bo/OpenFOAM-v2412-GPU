/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::faceZoneToCell

Description
    A \c topoSetCellSource to select cells based on \c front
    or \c back side of given \c faceZone(s).

    Operands:
    \table
      Operand   | Type        | Location
      input     | faceZone(s) | $FOAM_CASE/constant/polyMesh/faceZones
      output    | cellSet     | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      faceZoneToCell;
        option      <option>;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        zones
        (
            <faceZoneName1>
            <faceZoneName2>
            ...
        );

        // Option-2
        zone <faceZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: faceZoneToCell         | word |  yes  | -
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
      front   | Cells on the front side of the faceZone (also as 'master')
      back    | Cells on the back side of the faceZone (also as 'slave')
      both    | Both front and back side of the faceZone
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                    | Type     | Req'd  | Dflt
      zones    | Names of input faceZones       | wordRes  | cond'l | -
      zone     | Name of input faceZone         | wordRe   | cond'l | -
    \endverbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c zones, and \c zone.

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource
    - Foam::faceToCell

SourceFiles
    faceZoneToCell.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneToCell_H
#define faceZoneToCell_H

#include "topoSetCellSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceZoneToCell Declaration
\*---------------------------------------------------------------------------*/

class faceZoneToCell
:
    public topoSetCellSource
{
public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            FRONT,  //!< The front (positive normal) side of the faces
            BACK,   //!< The back (negative normal) side of the faces
            BOTH,   //!< Both front and back side of faces
            // Compatibility
            MASTER = FRONT, SLAVE = BACK
        };


private:

    // Private Data

        static const Enum<faceAction> faceActionNames_;

        //- Add usage string
        static addToUsageTable usage_;

        //- Matcher for face zones
        wordRes zoneMatcher_;

        //- Selection type
        faceAction option_;


    // Private Member Functions

        void combine
        (
            topoSet& set,
            const labelUList& zoneIDs,
            const bool add,
            const bool verbosity
        ) const;

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("faceZoneToCell");


    // Constructors

        //- Construct from mesh, zones selector and selection option
        faceZoneToCell
        (
            const polyMesh& mesh,
            const wordRes& zoneSelector,
            const faceAction option
        );

        //- Construct from mesh, single zone selector and selection option
        faceZoneToCell
        (
            const polyMesh& mesh,
            const wordRe& zoneName,
            const faceAction option
        );

        //- Construct from dictionary
        faceZoneToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceZoneToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceZoneToCell() = default;


    // Member Functions

        //- Return the current zones selector
        const wordRes& zones() const noexcept;

        //- Define the zones selector
        void zones(const wordRes& zoneSelector);

        //- Define the zones selector with a single zone selector
        void zones(const wordRe& zoneName);


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
