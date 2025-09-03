/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::zoneToFace

Description
    A \c topoSetFaceSource to convert \c faceZone(s) to a \c faceSet.

    Operands:
    \table
      Operand   | Type        | Location
      input     | faceZone(s) | $FOAM_CASE/constant/polyMesh/faceZones
      output    | faceSet     | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      zoneToFace;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        zones
        (
            <faceZoneName0>
            <faceZoneName1>
            ...
        );

        // Option-2
        zone     <faceZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: zoneToFace             | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      zones    | Names of input faceZones   | wordRes  | cond'l | -
      zone     | Name of input faceZone     | wordRe   | cond'l | -
    \endverbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c zones, and \c zone.

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    zoneToFace.C

\*---------------------------------------------------------------------------*/

#ifndef zoneToFace_H
#define zoneToFace_H

#include "topoSetFaceSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneToFace Declaration
\*---------------------------------------------------------------------------*/

class zoneToFace
:
    public topoSetFaceSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Matcher for zones
        wordRes zoneMatcher_;

        //- Explicitly specified zone ids
        labelList zoneIDs_;



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
    TypeName("zoneToFace");


    // Constructors

        //- Construct from mesh and zones selector
        zoneToFace(const polyMesh& mesh, const wordRes& zoneSelector);

        //- Construct from mesh and single zone selector
        zoneToFace(const polyMesh& mesh, const wordRe& zoneName);

        //- Construct from mesh and specified zone IDs
        zoneToFace(const polyMesh& mesh, const labelUList& zoneIDs);

        //- Construct from dictionary
        zoneToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        zoneToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~zoneToFace() = default;


    // Member Functions

        //- Return the current zones selector
        const wordRes& zones() const noexcept;

        //- Define the zones selector
        void zones(const wordRes& zoneSelector);

        //- Define the zones selector with a single zone selector
        void zones(const wordRe& zoneName);

        //- Define the faceZone IDs to use (must exist).
        //  Clears the zone name matcher
        void zones(const labelUList& zoneIDs);

        //- Define the faceZone ID to use (must exist).
        //  Clears the zone name matcher
        void zones(const label zoneID);


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
