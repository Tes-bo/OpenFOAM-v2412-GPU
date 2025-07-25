/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::surfaceToPoint

Description
    A \c topoSetPointSource to select points based on
    relation to a surface given by an external file.

    \c surfaceToPoint can select based on:
    - distance to surface
    - inside/outside status to surface
      (Uses normal of nearest surface triangle so requires valid surface
      topology.)

    Operands:
    \table
      Operand    | Type       | Location
      input      | triSurface | $FOAM_CASE/constant/triSurface/\<file\>
      output     | pointSet   | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name            <name>;
        type            pointSet;
        action          <action>;

        // Mandatory entries
        source          surfaceToPoint;
        file            <surfaceFileName>;
        includeInside   false;
        includeOutside  true;
        nearDistance    0.5;

        // Optional entries
        fileType        stl;
        scale           2.0;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                           | Type | Req'd | Dflt
      name       | Name of pointSet                      | word |  yes  | -
      type       | Type name: pointSet                   | word |  yes  | -
      action     | Action applied on points - see below  | word |  yes  | -
      source     | Source name: surfaceToPoint           | word |  yes  | -
      file       | The surface "filename"                | word |  yes  | -
      fileType   | The format of the surface file        | word |  no   | ""
      nearDistance | Near distance to the surface        | scalar | yes | -
      includeInside | Flag to include inside cells       | bool |  yes  | -
      includeOutside | Flag to include outside cells     | bool |  yes  | -
      scale        | Surface scaling factor              | scalar | no  | -1
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointSet from selected points
      add      | Add selected points into this pointSet
      subtract | Remove selected points from this pointSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    surfaceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceToPoint_H
#define surfaceToPoint_H

#include "topoSetPointSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceToPoint Declaration
\*---------------------------------------------------------------------------*/

class surfaceToPoint
:
    public topoSetPointSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of surface file
        const fileName surfName_;

        //- Surface file type (optional)
        const word surfType_;

        //- Surface scaling factor (optional)
        const scalar scale_;

        //- If > 0 : include points with distance to surface less than nearDist.
        const scalar nearDist_;

        //- Include all points on opposite sign of normal on nearest surface.
        const bool includeInside_;

        //- Include all points on this sign of normal on nearest surface.
        const bool includeOutside_;


    // Private Member Functions

        //- Depending on surface add to or delete from pointSet.
        void combine(topoSet& set, const bool add) const;

        //- Check settings at construction time.
        void checkSettings() const;


public:

    //- Runtime type information
    TypeName("surfaceToPoint");


    // Constructors

        //- Construct from components
        surfaceToPoint
        (
            const polyMesh& mesh,
            const fileName& surfName,
            const scalar nearDist,
            const bool includeInside,
            const bool includeOutside
        );

        //- Construct from dictionary
        surfaceToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        surfaceToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~surfaceToPoint() = default;


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
