/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::faceBoolSet

Description
    A special purpose topoSet with the face labels stored as a boolList.
    It does not correspond to a faceSet either (no associated IOobject).

SourceFiles
    faceBoolSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceBoolSet_H
#define Foam_faceBoolSet_H

#include "topoBoolSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceBoolSet Declaration
\*---------------------------------------------------------------------------*/

class faceBoolSet
:
    public topoBoolSet
{
public:

    //- Runtime type information
    TypeNameNoDebug("faceBoolSet");


    // Constructors

        //- Construct with nFaces elements, all elements unset
        explicit faceBoolSet(const polyMesh& mesh);

        //- Construct with nFaces elements, using initial val
        faceBoolSet(const polyMesh& mesh, const bool val);

        //- Copy construct from list, resizing to nFaces elements as required
        faceBoolSet(const polyMesh& mesh, const boolList& bools);

        //- Move construct from list, resizing to nFaces elements as required
        faceBoolSet(const polyMesh& mesh, boolList&& bools);


    //- Destructor
    virtual ~faceBoolSet() = default;


    // Member Functions

        //- Sync faceBoolSet across coupled patches.
        virtual void sync(const polyMesh& mesh);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels.
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const mapDistributePolyMesh& map);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh& mesh,
            const label maxLen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
