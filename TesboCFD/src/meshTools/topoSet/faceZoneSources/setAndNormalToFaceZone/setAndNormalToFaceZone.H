/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::setAndNormalToFaceZone

Description
    A \c topoSetSource to select faces based on usage in a \c faceSet,
    where the normal vector is used to orient the faces.

    Operands:
    \table
      Operand   | Type     | Location
      input     | faceSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      setAndNormalToFaceZone;
        faceSet     <faceSetName>;
        normal      (1 0 0);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceZoneSet                 | word |  yes  | -
      type       | Type name: faceZoneSet              | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: setAndNormalToFaceZone | word |  yes  | -
      faceSet    | Name of input faceSet               | word |  yes  | -
      normal     | The normal vector for orientation   | vector | yes | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceZoneSet from selected faces
      add      | Add selected faces into this faceZoneSet
      subtract | Remove selected faces from this faceZoneSet
    \endverbatim

SourceFiles
    setAndNormalToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setAndNormalToFaceZone_H
#define setAndNormalToFaceZone_H

#include "topoSetFaceZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class setAndNormalToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setAndNormalToFaceZone
:
    public topoSetFaceZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;

        //- Normal used to orient the faces
        vector normal_;


public:

    //- Runtime type information
    TypeName("setAndNormalToFaceZone");

    // Constructors

        //- Construct from components
        setAndNormalToFaceZone
        (
            const polyMesh& mesh,
            const word& setName,
            const vector& normal
        );

        //- Construct from dictionary
        setAndNormalToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setAndNormalToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setAndNormalToFaceZone() = default;


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
