/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::hexRef8Data

Description
    Various for reading/decomposing/reconstructing/distributing refinement
    data.

SourceFiles
    hexRef8Data.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_hexRef8Data_H
#define Foam_hexRef8Data_H

#include "labelIOList.H"
#include "uniformDimensionedFields.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class mapDistributePolyMesh;
class refinementHistory;
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class hexRef8Data Declaration
\*---------------------------------------------------------------------------*/

class hexRef8Data
{
    // Private Data

        autoPtr<labelIOList> cellLevelPtr_;

        autoPtr<labelIOList> pointLevelPtr_;

        autoPtr<uniformDimensionedScalarField> level0EdgePtr_;

        autoPtr<refinementHistory> refHistoryPtr_;


public:

    // Generated Methods

        //- No copy construct
        hexRef8Data(const hexRef8Data&) = delete;

        //- No copy assignment
        void operator=(const hexRef8Data&) = delete;


    // Constructors

        //- Construct from IOobject, reading if MUST_READ or LAZY_READ.
        //- Use sync() if used in redistributing to different number of ranks
        explicit hexRef8Data(const IOobject& io);

        //- Construct as subset
        hexRef8Data
        (
            const IOobject& io,
            const hexRef8Data&,
            const labelList& cellMap,
            const labelList& pointMap
        );

        //- Construct from multiple hexRef8Data
        hexRef8Data
        (
            const IOobject& io,
            const UPtrList<const labelList>& cellMaps,
            const UPtrList<const labelList>& pointMaps,
            const UPtrList<const hexRef8Data>&
        );


    //- Destructor
    ~hexRef8Data();


    // Member Functions

        //- Parallel synchronise. This enforces valid objects on all processors
        //- (even if they don't have a mesh). Used by redistributePar.
        void sync(const IOobject& io);

        //- Update local numbering for changed mesh.
        void updateMesh(const mapPolyMesh&);

        //- In-place distribute
        void distribute(const mapDistributePolyMesh&);

        //- Write
        bool write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
