/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::oversetGAMGInterfaceField

Description
    GAMG agglomerated

SourceFiles
    oversetGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef oversetGAMGInterfaceField_H
#define oversetGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "oversetGAMGInterface.H"
#include "oversetLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class oversetGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class oversetGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public oversetLduInterfaceField
{
    // Private Member Functions

        //- No copy construct
        oversetGAMGInterfaceField
        (
            const oversetGAMGInterfaceField&
        ) = delete;

        //- No copy assignment
        void operator=(const oversetGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        oversetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        oversetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );

        //- Construct from GAMG interface and Istream
        oversetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            Istream& is
        );

        //- Construct by assembling and return a clone.
        virtual autoPtr<GAMGInterfaceField> clone
        (
            const GAMGInterface& GAMGCp,
            const UPtrList<lduInterfaceField>& other
        ) const
        {
            NotImplemented;
            return nullptr;
        }


    //- Destructor
    virtual ~oversetGAMGInterfaceField();


    // Member Functions

        // Interface matrix update

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
