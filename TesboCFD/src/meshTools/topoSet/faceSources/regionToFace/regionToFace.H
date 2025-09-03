/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::regionToFace

Description
    A \c topoSetFaceSource to select cells belonging to a topologically
    connected region (that contains given points).

    Operands:
    \table
      Operand   | Type     | Location
      input     | region   | $FOAM_CASE/constant/\{\<region\>, polyMesh\}
      output    | faceSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      regionToFace;
        set         <faceSetName>;
        nearPoint   <point>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: regionToFace           | word |  yes  | -
      set / zone | Name of faceSet / zone restricting search  | word |  yes  | -
      nearPoint  | The point on/near to the region     | vector | yes | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    regionToFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionToFace_H
#define Foam_regionToFace_H

#include "topoSetFaceSource.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class regionToFace Declaration
\*---------------------------------------------------------------------------*/

class regionToFace
:
    public topoSetFaceSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set or zone to use
        word setName_;

        //- Is name a set or a zone
        bool isZone_;

        //- Coordinate that is nearest/on connected region
        const point nearPoint_;


    // Private Member Functions

        //- Walk edge-face-edge
        void markZone
        (
            const indirectPrimitivePatch& patch,
            const label proci,
            const label facei,
            const label zoneI,
            labelList& faceZone
        ) const;

        void combine
        (
            topoSet& set,
            const bool add,
            const labelUList& ids
        ) const;


public:

    //- Runtime type information
    TypeName("regionToFace");


    // Constructors

        //- Construct from components
        regionToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const point& nearPoint
        );

        //- Construct from dictionary
        regionToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        regionToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~regionToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
