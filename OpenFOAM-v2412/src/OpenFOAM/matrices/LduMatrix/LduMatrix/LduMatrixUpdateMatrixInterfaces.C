/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "LduMatrix.H"
#include "lduInterfaceField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::initMatrixInterfaces
(
    const bool add,
    const FieldField<Field, LUType>& interfaceCoeffs,
    const Field<Type>& psiif,
    Field<Type>& result
) const
{
    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces_, interfacei)
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    lduMesh_.lduAddr(),
                    interfacei,
                    psiif,
                    interfaceCoeffs[interfacei],
                    commsType
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces_.size();
            interfacei++
        )
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    lduMesh_.lduAddr(),
                    interfacei,
                    psiif,
                    interfaceCoeffs[interfacei],
                    UPstream::commsTypes::buffered
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << exit(FatalError);
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, LUType>& interfaceCoeffs,
    const Field<Type>& psiif,
    Field<Type>& result,
    const label startRequest
) const
{
    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    if
    (
        commsType == UPstream::commsTypes::nonBlocking
     && UPstream::nPollProcInterfaces
    )
    {
        // Wait for some interface requests to become available and
        // consume them. No guarantee that the finished requests actually
        // correspond to any particular interface, but it is reasonably
        // probable that some interfaces will be able to start consumption
        // without waiting for all requests.

        DynamicList<int> indices;  // (work array)

        for
        (
            bool pollingActive = (UPstream::nPollProcInterfaces < 0);
            (
                pollingActive
             && UPstream::waitSomeRequests(startRequest, -1, &indices)
            );
            /*nil*/
        )
        {
            pollingActive = false;

            forAll(interfaces_, interfacei)
            {
                auto* intf = interfaces_.get(interfacei);

                if (intf && !intf->updatedMatrix())
                {
                    if (intf->ready())
                    {
                        intf->updateInterfaceMatrix
                        (
                            result,
                            add,
                            lduMesh_.lduAddr(),
                            interfacei,
                            psiif,
                            interfaceCoeffs[interfacei],
                            commsType
                        );
                    }
                    else
                    {
                        pollingActive = true;
                    }
                }
            }
        }
    }


    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        // Wait until sends/receives have finished.
        // - effectively a no-op (without waiting) if already completed.
        if (commsType == UPstream::commsTypes::nonBlocking)
        {
            UPstream::waitRequests(startRequest);
        }

        // Check/no-check for updatedMatrix() ?
        const bool noCheck = (commsType == UPstream::commsTypes::buffered);

        forAll(interfaces_, interfacei)
        {
            auto* intf = interfaces_.get(interfacei);

            if (intf && (noCheck || !intf->updatedMatrix()))
            {
                intf->updateInterfaceMatrix
                (
                    result,
                    add,
                    lduMesh_.lduAddr(),
                    interfacei,
                    psiif,
                    interfaceCoeffs[interfacei],
                    commsType
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        for (const auto& schedEval : patchSchedule)
        {
            const label interfacei = schedEval.patch;

            if (interfaces_.set(interfacei))
            {
                if (schedEval.init)
                {
                    interfaces_[interfacei].initInterfaceMatrixUpdate
                    (
                        result,
                        add,
                        lduMesh_.lduAddr(),
                        interfacei,
                        psiif,
                        interfaceCoeffs[interfacei],
                        commsType
                    );
                }
                else
                {
                    interfaces_[interfacei].updateInterfaceMatrix
                    (
                        result,
                        add,
                        lduMesh_.lduAddr(),
                        interfacei,
                        psiif,
                        interfaceCoeffs[interfacei],
                        commsType
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces_.size();
            interfacei++
        )
        {
            if (interfaces_.set(interfacei))
            {
                interfaces_[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    lduMesh_.lduAddr(),
                    interfacei,
                    psiif,
                    interfaceCoeffs[interfacei],
                    UPstream::commsTypes::buffered
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
