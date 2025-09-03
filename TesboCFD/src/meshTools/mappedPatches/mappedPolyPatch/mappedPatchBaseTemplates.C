/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

template<class Type>
void Foam::mappedPatchBase::distribute(List<Type>& lst) const
{
    const label myComm = getCommunicator();  // Get or create

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const auto& interp = AMI();
            const label oldWarnComm = UPstream::commWarn(myComm);
            const label oldWorldComm = UPstream::commWorld(myComm);

            if (sameWorld())
            {
                // lst is the other side's values
                lst = interp.interpolateToSource(Field<Type>(std::move(lst)));
            }
            else
            {
                // lst is my local data. Now the mapping in the AMI is
                // from my side to other side. Each processor contains either
                // faces from one side or from the other side.

                if (masterWorld())
                {
                    // I have lst.size() faces on my side, zero of the other
                    // side

                    tmp<Field<Type>> tmasterFld
                    (
                        interp.interpolateToSource(Field<Type>(0))
                    );
                    (void)interp.interpolateToTarget
                    (
                        Field<Type>(std::move(lst))
                    );

                    // We've received in our interpolateToSource the
                    // contribution from the other side
                    lst = tmasterFld;
                }
                else
                {
                    (void)interp.interpolateToSource
                    (
                        Field<Type>(std::move(lst))
                    );
                    tmp<Field<Type>> tmasterFld
                    (
                        interp.interpolateToTarget(Field<Type>(0))
                    );

                    // We've received in our interpolateToTarget the
                    // contribution from the other side
                    lst = tmasterFld;
                }
            }

            // Restore communicator settings
            UPstream::commWarn(oldWarnComm);
            UPstream::commWorld(oldWorldComm);
            break;
        }
        default:
        {
            const auto& m = map();
            const label oldWarnComm = UPstream::commWarn(m.comm());

            m.distribute(lst);

            UPstream::commWarn(oldWarnComm);
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::distribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    const label myComm = getCommunicator();  // Get or create

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const auto& interp = AMI();
            const label oldWarnComm = UPstream::commWarn(myComm);
            const label oldWorldComm = UPstream::commWorld(myComm);

            lst = interp.interpolateToSource(Field<Type>(std::move(lst)), cop);

            UPstream::commWarn(oldWarnComm);
            UPstream::commWorld(oldWorldComm);
            break;
        }
        default:
        {
            // Force early construction of parallel data
            (void)patch_.boundaryMesh().mesh().tetBasePtIs();
            const auto& m = map();

            const label oldWarnComm = UPstream::commWarn(myComm);
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                m.schedule(),
                m.constructSize(),
                m.subMap(),
                false,
                m.constructMap(),
                false,
                lst,
                Type(Zero),
                cop,
                flipOp(),
                UPstream::msgType(),
                myComm
            );
            UPstream::commWarn(oldWarnComm);
        }
    }
}


template<class Type>
void Foam::mappedPatchBase::reverseDistribute(List<Type>& lst) const
{
    const label myComm = getCommunicator();  // Get or create

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const auto& interp = AMI();
            const label oldWarnComm = UPstream::commWarn(myComm);
            const label oldWorldComm = UPstream::commWorld(myComm);

            lst = interp.interpolateToTarget(Field<Type>(std::move(lst)));

            UPstream::commWarn(oldWarnComm);
            UPstream::commWorld(oldWorldComm);
            break;
        }
        default:
        {
            // Force early construction of parallel data
            (void)patch_.boundaryMesh().mesh().tetBasePtIs();
            const auto& m = map();

            const label oldWarnComm = UPstream::commWarn(m.comm());
            m.reverseDistribute(sampleSize(), lst);
            UPstream::commWarn(oldWarnComm);
            break;
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::reverseDistribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    const label myComm = getCommunicator();  // Get or create

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const auto& interp = AMI();
            const label oldWarnComm = UPstream::commWarn(myComm);
            const label oldWorldComm = UPstream::commWorld(myComm);

            lst = interp.interpolateToTarget(Field<Type>(std::move(lst)), cop);

            UPstream::commWarn(oldWarnComm);
            UPstream::commWorld(oldWorldComm);
            break;
        }
        default:
        {
            (void)patch_.boundaryMesh().mesh().tetBasePtIs();
            const auto& m = map();
            const label cSize = sampleSize();

            const label oldWarnComm = UPstream::commWarn(myComm);

            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                m.schedule(),
                cSize,
                m.constructMap(),
                false,
                m.subMap(),
                false,
                lst,
                Type(Zero),
                cop,
                flipOp(),
                UPstream::msgType(),
                myComm
            );

            UPstream::commWarn(oldWarnComm);
            break;
        }
    }
}


template<class Type>
bool Foam::mappedPatchBase::writeIOField
(
    const regIOobject& obj,
    dictionary& dict
)
{
    const auto* fldPtr = isA<IOField<Type>>(obj);
    if (fldPtr)
    {
        const auto& fld = *fldPtr;

        primitiveEntry* pePtr = new primitiveEntry
        (
            fld.name(),
            token(new token::Compound<List<Type>>(fld))
        );

        dict.set(pePtr);
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
bool Foam::mappedPatchBase::constructIOField
(
    const word& name,
    token& tok,
    Istream& is,
    objectRegistry& obr
)
{
    if (tok.isCompound<List<Type>>())
    {
        auto* fldPtr = obr.getObjectPtr<IOField<Type>>(name);
        if (!fldPtr)
        {
            fldPtr = new IOField<Type>
            (
                IOobject
                (
                    name,
                    obr,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::REGISTER
                ),
                Foam::zero{}
            );
            regIOobject::store(fldPtr);
        }

        fldPtr->transfer
        (
            tok.transferCompoundToken<List<Type>>(is)
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::mappedPatchBase::storeField
(
    objectRegistry& obr,
    const word& fieldName,
    const Field<Type>& values
)
{
    auto* fldPtr = obr.getObjectPtr<IOField<Type>>(fieldName);
    if (!fldPtr)
    {
        fldPtr = new IOField<Type>
        (
            IOobject
            (
                fieldName,
                obr,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            Foam::zero{}
        );
        regIOobject::store(fldPtr);
    }

    *fldPtr = values;
}


// ************************************************************************* //
