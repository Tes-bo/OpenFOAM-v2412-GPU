/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::setsToFaceZone

Description
    A \c topoSetSource to select all faces in a given \c faceSet,
    which are oriented in slave cells of a given \c cellSet.

    Operands:
    \table
      Operand   | Type     | Location
      input 1   | faceSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      input 2   | cellSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | faceZone | $FOAM_CASE/constant/polyMesh/faceZones
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceZoneSet;
        action      <action>;

        // Mandatory entries
        source      setToFaceZone;
        faceSet     <faceSetName>;
        cellSet     <cellSetName>;

        // Optional entries
        flip        false;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceZone                    | word |  yes  | -
      type       | Type name: faceZoneSet              | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: setsToFaceZone         | word |  yes  | -
      faceSet    | Name of input faceSet               | word |  yes  | -
      cellSet    | Name of input cellSet containing the slave cells <!--
                                                   --> | word |  yes  | -
      flip       | Flag to select master/slave cells   | bool |  no   | false
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceZone from selected faces
      add      | Add selected faces of a faceZoneSet into this faceZone
      subtract | Remove selected faces of a faceZoneSet from this faceZone
    \endverbatim

See also
    - Foam::topoSetSource

SourceFiles
    setsToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setsToFaceZone_H
#define setsToFaceZone_H

#include "topoSetFaceZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setsToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setsToFaceZone
:
    public topoSetFaceZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        const word faceSetName_;

        //- Name of set to use
        const word cellSetName_;

        //- Whether cellSet is slave cells or master cells
        const bool flip_;


public:

    //- Runtime type information
    TypeName("setsToFaceZone");


    // Constructors

        //- Construct from components
        setsToFaceZone
        (
            const polyMesh& mesh,
            const word& faceSetName,
            const word& cellSetName,
            const bool flip
        );

        //- Construct from dictionary
        setsToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setsToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setsToFaceZone() = default;


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
