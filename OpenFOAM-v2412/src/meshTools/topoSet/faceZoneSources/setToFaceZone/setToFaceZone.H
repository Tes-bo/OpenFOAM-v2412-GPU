/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::setToFaceZone

Description
    A \c topoSetSource to convert a \c faceSet
    to a \c faceZone (and associated \c faceSet).

    Operands:
    \table
      Operand   | Type     | Location
      input     | faceSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 1  | faceSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 2  | faceZone | $FOAM_CASE/constant/polyMesh/faceZones
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
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceZone                    | word |  yes  | -
      type       | Type name: faceZoneSet              | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: setToFaceZone          | word |  yes  | -
      faceSet    | Name of input faceSet               | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceZone from selected faces of a faceSet
      add      | Add selected faces of a faceSet into this faceZone
      subtract | Remove selected faces of a faceSet from this faceZone
    \endverbatim

Note
    The use of keyword \c faceSet is consistent with the \c setsToFaceZone,
    but inconsistent with \c setToCellZone and \c setToPointZone.
    Both of which use \c set instead.

See also
    - Foam::topoSetSource

SourceFiles
    setToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setToFaceZone_H
#define setToFaceZone_H

#include "topoSetFaceZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setToFaceZone
:
    public topoSetFaceZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;


public:

    //- Runtime type information
    TypeName("setToFaceZone");


    // Constructors

        //- Construct from components
        setToFaceZone
        (
            const polyMesh& mesh,
            const word& setName
        );

        //- Construct from dictionary
        setToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setToFaceZone() = default;


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
