/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "lduPrimitiveProcessorInterface.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(lduPrimitiveProcessorInterface);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveProcessorInterface::lduPrimitiveProcessorInterface
(
    const lduPrimitiveProcessorInterface& p
)
:
    faceCells_(p.faceCells()),
    myProcNo_(p.myProcNo()),
    neighbProcNo_(p.neighbProcNo()),
    forwardT_(p.forwardT()),
    tag_(p.tag()),
    comm_(p.comm())
{}


Foam::lduPrimitiveProcessorInterface::lduPrimitiveProcessorInterface
(
    const labelUList& faceCells,
    const label myProcNo,
    const label neighbProcNo,
    const tensorField& forwardT,
    const int tag,
    const label comm
)
:
    faceCells_(faceCells),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    forwardT_(forwardT),
    tag_(tag),
    comm_(comm)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField>
Foam::lduPrimitiveProcessorInterface::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return interfaceInternalField(internalData, faceCells_);
}


Foam::tmp<Foam::labelField>
Foam::lduPrimitiveProcessorInterface::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    auto tfld = tmp<labelField>::New(faceCells.size());
    auto& fld = tfld.ref();

    forAll(faceCells, i)
    {
        fld[i] = internalData[faceCells[i]];
    }
    return tfld;
}


void Foam::lduPrimitiveProcessorInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    processorLduInterface::send(commsType, interfaceInternalField(iF)());
}


void Foam::lduPrimitiveProcessorInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF,
    const labelUList& faceCells
) const
{
    processorLduInterface::send
    (
        commsType,
        interfaceInternalField(iF, faceCells)()
    );
}


Foam::tmp<Foam::labelField>
Foam::lduPrimitiveProcessorInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return processorLduInterface::receive<label>(commsType, faceCells_.size());
}


// ************************************************************************* //
