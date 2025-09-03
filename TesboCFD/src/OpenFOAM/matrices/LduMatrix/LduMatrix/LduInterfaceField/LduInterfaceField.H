/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::LduInterfaceField

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

SourceFiles
    LduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LduInterfaceField_H
#define Foam_LduInterfaceField_H

#include "lduInterfaceField.H"
#include "primitiveFieldsFwd.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class lduMatrix;
template<class Type> class LduInterfaceField;

/*---------------------------------------------------------------------------*\
                     Class LduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LduInterfaceField
:
    public lduInterfaceField
{
    // Private Member Functions

        //- No copy construct
        LduInterfaceField(const LduInterfaceField&) = delete;

        //- No copy assignment
        void operator=(const LduInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("LduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        LduInterfaceField(const lduInterface& patch)
        :
            lduInterfaceField(patch)
        {}


    //- Destructor
    virtual ~LduInterfaceField() = default;


    // Member Functions

        // Coupled interface functionality

            //- Inherit initInterfaceMatrixUpdate from lduInterfaceField
            using lduInterfaceField::initInterfaceMatrixUpdate;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>&,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Inherit updateInterfaceMatrix from lduInterfaceField
            using lduInterfaceField::updateInterfaceMatrix;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
