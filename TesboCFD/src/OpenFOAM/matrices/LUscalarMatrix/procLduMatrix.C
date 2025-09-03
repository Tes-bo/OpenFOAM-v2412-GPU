/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "procLduMatrix.H"
#include "procLduInterface.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::procLduMatrix::procLduMatrix
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    upperAddr_(ldum.lduAddr().upperAddr()),
    lowerAddr_(ldum.lduAddr().lowerAddr()),
    diag_(ldum.diag()),
    upper_(ldum.upper()),
    lower_(ldum.lower()),
    interfaces_(interfaces.count_nonnull())
{
    label nInterfaces = 0;

    forAll(interfaces, i)
    {
        if (interfaces.test(i))
        {
            interfaces_.set
            (
                nInterfaces++,
                new procLduInterface
                (
                    interfaces[i],
                    interfaceCoeffs[i]
                )
            );
        }
    }
}


Foam::procLduMatrix::procLduMatrix(Istream& is)
{
    read(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::procLduMatrix::clear()
// {
//     upperAddr_.clear();
//     lowerAddr_.clear();
//     diag_.clear();
//     upper_.clear();
//     lower_.clear();
//     interfaces_.clear();
// }


void Foam::procLduMatrix::read(Istream& is)
{
    is >> upperAddr_ >> lowerAddr_ >> diag_ >> upper_ >> lower_ >> interfaces_;
}


void Foam::procLduMatrix::write(Ostream& os) const
{
    os << upperAddr_ << lowerAddr_ << diag_ << upper_ << lower_ << interfaces_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, procLduMatrix& mat)
{
    mat.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const procLduMatrix& mat)
{
    mat.write(os);
    return os;
}


// ************************************************************************* //
