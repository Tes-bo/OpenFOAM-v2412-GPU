/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::lduInterface

Description
    An abstract base class for implicitly-coupled interfaces
    e.g. processor and cyclic patches.

SourceFiles
    lduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduInterface_H
#define Foam_lduInterface_H

#include "labelField.H"
#include "typeInfo.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class lduInterface;

/*---------------------------------------------------------------------------*\
                        Class lduInterface Declaration
\*---------------------------------------------------------------------------*/

class lduInterface
{
public:

    //- Runtime type information
    TypeNameNoDebug("lduInterface");


    // Constructors

        //- Default construct
        lduInterface() noexcept = default;


    //- Destructor
    virtual ~lduInterface() = default;


    // Member Functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const = 0;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const = 0;

            //- Return the values of the given internal data adjacent to
            //- the interface as a field using faceCell mapping
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData,
                const labelUList& faceCells
            ) const = 0;

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const
            {}

            //- Initialise transfer of internal field adjacent to the interface
            //- using faceCells mapping
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF,
                const labelUList& faceCells
            ) const
            {}

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
