/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::eagerGAMGProcAgglomeration

Description
    'Eager' processor agglomeration of GAMGAgglomerations: at every
    level agglomerates 'mergeLevels' number of processors onto the
    minimum processor number.

SourceFiles
    eagerGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_eagerGAMGProcAgglomeration_H
#define Foam_eagerGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class eagerGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class eagerGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- Agglpmeration level
        const label mergeLevels_;


public:

    //- Runtime type information
    TypeName("eager");

    //- No copy construct
    eagerGAMGProcAgglomeration(const eagerGAMGProcAgglomeration&) = delete;

    //- No copy assignment
    void operator=(const eagerGAMGProcAgglomeration&) = delete;


    // Constructors

        //- Construct given agglomerator and controls
        eagerGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~eagerGAMGProcAgglomeration();


    // Member Functions

        //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
