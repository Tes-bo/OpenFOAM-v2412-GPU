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
    Foam::faceToFace

Description
    A \c topoSetFaceSource to select faces based on usage in another \c faceSet.

    Operands:
    \table
      Operand   | Type    | Location
      input     | faceSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | faceSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      faceToFace;

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
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: faceToFace             | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      sets     | Names of input faceSets    | wordList | cond'l | -
      zones    | Names of input faceZones   | wordList | cond'l | -
      set      | Name of input faceSet      | word     | cond'l | -
      zone     | Name of input faceZone     | word     | cond'l | -
    \endverbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c sets, and \c set.

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    faceToFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceToFace_H
#define Foam_faceToFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceToFace Declaration
\*---------------------------------------------------------------------------*/

class faceToFace
:
    public topoSetFaceSource
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
    TypeName("faceToFace");


    // Constructors

        //- Construct from components
        faceToFace(const polyMesh& mesh, const word& setName);

        //- Construct from dictionary
        faceToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceToFace() = default;


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
