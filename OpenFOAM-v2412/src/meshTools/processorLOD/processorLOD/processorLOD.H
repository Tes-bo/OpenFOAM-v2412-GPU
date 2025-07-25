/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::processorLOD

Description
    Base class to generate a parallel distribution map for sending sufficient
    target objects to cover a description of the source object, based on
    processor Level Of Detail (LOD) shapes

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorLOD_H
#define Foam_processorLOD_H

#include "autoPtr.H"
#include "typeInfo.H"
#include "mapDistributeBase.H"  // For layoutTypes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistribute;

/*---------------------------------------------------------------------------*\
                        Class processorLOD Declaration
\*---------------------------------------------------------------------------*/

class processorLOD
{
protected:

    // Protected Data

        //- Maximum number of objects per leaf
        label maxObjectsPerLeaf_;

        //- Number of objects of this type.
        //- e.g. number of faces/cells on this processor
        label nObjectsOfType_;


public:

    //- Runtime type information
    TypeName("processorLOD");

    //- Construct from components
    processorLOD
    (
        const label maxObjectsPerLeaf,
        const label nObjectsOfType
    ) noexcept;

    //- Destructor
    virtual ~processorLOD() = default;


    // Member Functions

        //- Return the parallel distribution map
        //- (usually linear construct order)
        virtual autoPtr<mapDistribute> map
        (
            const mapDistributeBase::layoutTypes constructLayout
          = mapDistributeBase::layoutTypes::linear
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
