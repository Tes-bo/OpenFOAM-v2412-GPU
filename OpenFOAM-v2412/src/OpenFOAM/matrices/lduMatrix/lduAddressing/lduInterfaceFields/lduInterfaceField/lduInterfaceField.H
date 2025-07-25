/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::lduInterfaceField

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

SourceFiles
    lduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduInterfaceField_H
#define Foam_lduInterfaceField_H

#include "lduInterface.H"
#include "lduAddressing.H"
#include "primitiveFieldsFwd.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class lduMatrix;
class lduInterfaceField;

/*---------------------------------------------------------------------------*\
                     Class lduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class lduInterfaceField
{
    // Private Data

        //- Reference to the coupled patch this field is defined for
        const lduInterface& interface_;

        //- Update index used so that updateInterfaceMatrix is called only once
        //- during the construction of the matrix
        mutable bool updatedMatrix_;


public:

    //- Runtime type information
    TypeName("lduInterfaceField");


    // Generated Methods

        //- No copy construct
        lduInterfaceField(const lduInterfaceField&) = delete;

        //- No copy assignment
        void operator=(const lduInterfaceField&) = delete;


    // Constructors

        //- Construct given coupled patch
        explicit lduInterfaceField(const lduInterface& patch)
        :
            interface_(patch),
            updatedMatrix_(false)
        {}


    //- Destructor
    virtual ~lduInterfaceField() = default;


    // Member Functions

        //- Return the interface
        const lduInterface& interface() const noexcept
        {
            return interface_;
        }

        //- Return the interface type
        virtual const word& interfaceFieldType() const
        {
            return type();
        }


    // Coupled Interface

        //- Are all (receive) data available?
        virtual bool ready() const
        {
            return true;
        }

        //- Whether matrix has been updated
        bool updatedMatrix() const noexcept
        {
            return updatedMatrix_;
        }

        //- Set matrix as update-to-date, return the previous value
        bool updatedMatrix(bool flag) const noexcept
        {
            bool old(updatedMatrix_);
            updatedMatrix_ = flag;
            return old;
        }


        // Coupled interface matrix update

            //- Initialise neighbour matrix update.
            //- Add/subtract coupled contributions to matrix
            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality.
            //- Add/subtract coupled contributions to matrix
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const = 0;

            //- Add/subtract weighted contributions to internal field
            template<class Type>
            void addToInternalField
            (
                Field<Type>& result,
                const bool add,
                const labelUList& faceCells,
                const scalarField& coeffs,
                const Field<Type>& vals
            ) const;


    // Housekeeping

        //- Adjust whether matrix has been updated.
        //  \deprecated Prefer the updatedMatrix(bool) setter (JAN-2023)
        bool& updatedMatrix() noexcept
        {
            return updatedMatrix_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduInterfaceFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
