/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::procLduMatrix

Description
    I/O for lduMatrix and interface values.

SourceFiles
    procLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_procLduMatrix_H
#define Foam_procLduMatrix_H

#include "labelList.H"
#include "scalarField.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class procLduMatrix;
class procLduInterface;
class lduMatrix;

Istream& operator>>(Istream&, procLduMatrix&);
Ostream& operator<<(Ostream&, const procLduMatrix&);

/*---------------------------------------------------------------------------*\
                        Class procLduMatrix Declaration
\*---------------------------------------------------------------------------*/

class procLduMatrix
{
    // Private Data

        labelList upperAddr_;
        labelList lowerAddr_;
        scalarField diag_;
        scalarField upper_;
        scalarField lower_;
        PtrList<procLduInterface> interfaces_;

public:

    //- Friendship
    friend class LUscalarMatrix;


    // Generated Methods

        //- Default construct
        procLduMatrix() = default;

        //- No copy construct
        procLduMatrix(const procLduMatrix&) = delete;


    // Constructors

        //- Construct from components. Extracts active interfaces
        procLduMatrix
        (
            const lduMatrix& ldum,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Read construct from Istream
        explicit procLduMatrix(Istream& is);


    // Member Functions

        label size() const noexcept
        {
            return diag_.size();
        }

        // void clear();


    // IO Operations

        void read(Istream& is);
        void write(Ostream& os) const;

        friend Istream& operator>>(Istream&, procLduMatrix&);
        friend Ostream& operator<<(Ostream&, const procLduMatrix&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
