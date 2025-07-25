/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::clipPlaneToFace

Description
    A \c topoSetFaceSource to select all faces
    whose face centre is above the specified plane

    Operands:
    \table
      Operand   | Type    | Location
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
        source      clipPlaneToFace;

        point  (x y z);
        normal (x y z);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: clipPlaneToFace        | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                        | Type   | Req'd  | Dflt
      point    | The reference point for the plane  | vector | yes    | -
      normal   | The normal of the plane            | vector | yes    | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    clipPlaneToFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_clipPlaneToFace_H
#define Foam_clipPlaneToFace_H

#include "topoSetFaceSource.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class clipPlaneToFace Declaration
\*---------------------------------------------------------------------------*/

class clipPlaneToFace
:
    public topoSetFaceSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Point on the plane
        const vector point_;

        //- Normal to the plane
        const vector normal_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("clipPlaneToFace");


    // Constructors

        //- No default construct
        clipPlaneToFace() = delete;

        //- Construct from components
        clipPlaneToFace
        (
            const polyMesh& mesh,
            const point& basePoint,
            const vector& normal
        );

        //- Construct from dictionary
        clipPlaneToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from stream - a single box.
        clipPlaneToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~clipPlaneToFace() = default;


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
