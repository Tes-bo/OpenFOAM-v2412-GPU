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

#include "lduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduMatrix::initMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt
) const
{
    const UPstream::commsTypes commsType = UPstream::defaultCommsType;

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
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
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
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


void Foam::lduMatrix::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt,
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

        // const label maxPolling =
        // (
        //     (UPstream::nPollProcInterfaces < 0)
        //   ? (UPstream::nRequests() - startRequest)
        //   : UPstream::nPollProcInterfaces
        // );

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

            forAll(interfaces, interfacei)
            {
                auto* intf = interfaces.get(interfacei);

                if (intf && !intf->updatedMatrix())
                {
                    if (intf->ready())
                    {
                        intf->updateInterfaceMatrix
                        (
                            result,
                            add,
                            mesh().lduAddr(),
                            interfacei,
                            psiif,
                            coupleCoeffs[interfacei],
                            cmpt,
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

        // [OLDER]
        // Alternative for consuming interfaces as they become available.
        // Within the loop, the ready() calls an MPI_Test
        // that can trigger progression. However, a bit less reliably
        // (and less efficient) since it is implies multiple calls to
        // MPI_Test to progress the MPI transfer, but no guarantee that
        // any of them will actually complete within nPollProcInterfaces
        // calls.

        for (label i = 0; i < UPstream::nPollProcInterfaces; ++i)
        {
            bool allUpdated = true;

            forAll(interfaces, interfacei)
            {
                auto* intf = interfaces.get(interfacei);

                if (intf && !intf->updatedMatrix())
                {
                    if (intf->ready())
                    {
                        intf->updateInterfaceMatrix
                        (
                            result,
                            add,
                            mesh().lduAddr(),
                            interfacei,
                            psiif,
                            coupleCoeffs[interfacei],
                            cmpt,
                            commsType
                        );
                    }
                    else
                    {
                        allUpdated = false;
                    }
                }
            }

            if (allUpdated)
            {
                break;  // Early exit
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

        // Consume anything still outstanding
        forAll(interfaces, interfacei)
        {
            auto* intf = interfaces.get(interfacei);

            if (intf && (noCheck || !intf->updatedMatrix()))
            {
                intf->updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    commsType
                );
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        for (const auto& sched : patchSchedule)
        {
            const label interfacei = sched.patch;

            if (interfaces.set(interfacei))
            {
                if (sched.init)
                {
                    interfaces[interfacei].initInterfaceMatrixUpdate
                    (
                        result,
                        add,
                        mesh().lduAddr(),
                        interfacei,
                        psiif,
                        coupleCoeffs[interfacei],
                        cmpt,
                        commsType
                    );
                }
                else
                {
                    interfaces[interfacei].updateInterfaceMatrix
                    (
                        result,
                        add,
                        mesh().lduAddr(),
                        interfacei,
                        psiif,
                        coupleCoeffs[interfacei],
                        cmpt,
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
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
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
