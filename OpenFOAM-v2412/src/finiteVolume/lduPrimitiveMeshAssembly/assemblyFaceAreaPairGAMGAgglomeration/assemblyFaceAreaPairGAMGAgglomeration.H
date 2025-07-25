/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::assemblyFaceAreaPairGAMGAgglomeration

Description
    Agglomerate using the pair algorithm.

SourceFiles
    assemblyFaceAreaPairGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef assemblyFaceAreaPairGAMGAgglomeration_H
#define assemblyFaceAreaPairGAMGAgglomeration_H

#include "pairGAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class assemblyFaceAreaPairGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class assemblyFaceAreaPairGAMGAgglomeration
:
    public pairGAMGAgglomeration
{
public:

    //- Runtime type information
    TypeName("assembledFaceAreaPair");


    //- Destructor
    ~assemblyFaceAreaPairGAMGAgglomeration();


    // Constructors

        //- Construct given mesh and controls
        assemblyFaceAreaPairGAMGAgglomeration
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );

        //- Construct given mesh and controls
        assemblyFaceAreaPairGAMGAgglomeration
        (
            const lduMatrix& matrix,
            const scalarField& cellVolumes,
            const vectorField& faceAreas,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
