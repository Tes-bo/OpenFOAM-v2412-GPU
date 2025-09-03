/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "procLduInterface.H"
#include "lduInterfaceField.H"
#include "cyclicLduInterface.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::procLduInterface::procLduInterface
(
    const lduInterfaceField& interface,
    const scalarField& coeffs
)
:
    faceCells_(interface.interface().faceCells()),
    coeffs_(coeffs),
    myProcNo_(-1),
    neighbProcNo_(-1),
    tag_(-1),
    comm_(-1)
{
    if (isA<processorLduInterface>(interface.interface()))
    {
        const processorLduInterface& pldui =
            refCast<const processorLduInterface>(interface.interface());

        myProcNo_ = pldui.myProcNo();
        neighbProcNo_ = pldui.neighbProcNo();
        tag_ = pldui.tag();
        comm_ = pldui.comm();
    }
    else if (isA<cyclicLduInterface>(interface.interface()))
    {
    }
    else
    {
        FatalErrorInFunction
            << "Unknown lduInterface type "
            << interface.interface().type()
            << exit(FatalError);
    }
}


Foam::procLduInterface::procLduInterface(Istream& is)
{
    read(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::procLduInterface::read(Istream& is)
{
    is  >> faceCells_
        >> coeffs_
        >> myProcNo_
        >> neighbProcNo_
        >> tag_
        >> comm_;
}


void Foam::procLduInterface::write(Ostream& os) const
{
    os  << faceCells_
        << coeffs_
        << myProcNo_
        << neighbProcNo_
        << tag_
        << comm_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, procLduInterface& intf)
{
    intf.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const procLduInterface& intf)
{
    intf.write(os);
    return os;
}


// ************************************************************************* //
