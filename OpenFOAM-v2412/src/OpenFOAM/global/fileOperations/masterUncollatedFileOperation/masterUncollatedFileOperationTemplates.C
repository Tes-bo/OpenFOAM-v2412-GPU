/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "Pstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class FileOp>
Type Foam::fileOperations::masterUncollatedFileOperation::masterOp
(
    const fileName& fName,
    const FileOp& fop,
    const int tag,
    const label comm
) const
{
    if (IFstream::debug)
    {
        Pout<< "masterUncollatedFileOperation::masterOp : Operation "
            << typeid(FileOp).name()
            << " on " << fName << endl;
    }

    if (UPstream::is_parallel(comm))
    {
        const label myProci = UPstream::myProcNo(comm);
        const label numProc = UPstream::nProcs(comm);

        List<fileName> filePaths(numProc);
        filePaths[myProci] = fName;
        Pstream::gatherList(filePaths, tag, comm);
        // OR filePaths = Pstream::listGatherValues(fName, comm, tag)

        List<Type> result;
        if (UPstream::master(comm))
        {
            result.resize(numProc);
            result = fop(filePaths[0]);

            for (label i = 1; i < numProc; ++i)
            {
                if (filePaths[i] != filePaths[0])
                {
                    result[i] = fop(filePaths[i]);
                }
            }
        }

        return Pstream::listScatterValues(result, comm, tag);
    }

    return fop(fName);
}


template<class Type, class FileOp>
Type Foam::fileOperations::masterUncollatedFileOperation::masterOp
(
    const fileName& src,
    const fileName& dest,
    const FileOp& fop,
    const int tag,
    const label comm
) const
{
    if (IFstream::debug)
    {
        Pout<< "masterUncollatedFileOperation : Operation on src:" << src
            << " dest:" << dest << endl;
    }

    if (UPstream::is_parallel(comm))
    {
        const label myProci = UPstream::myProcNo(comm);
        const label numProc = UPstream::nProcs(comm);

        List<Pair<fileName>> filePaths(numProc);
        filePaths[myProci].first() = src;
        filePaths[myProci].second() = dest;
        Pstream::gatherList(filePaths, tag, comm);
        // OR
        // Pair<fileName> tup(src, dest);
        // filePaths = Pstream::listGatherValues(tup, comm, tag)

        List<Type> result;
        if (UPstream::master(comm))
        {
            result.resize(numProc);
            result = fop(filePaths[0].first(), filePaths[0].second());

            for (label i = 1; i < numProc; ++i)
            {
                // TBD: also check second() ?
                if (filePaths[i].first() != filePaths[0].first())
                {
                    result[i] =
                        fop(filePaths[i].first(), filePaths[i].second());
                }
            }
        }

        return Pstream::listScatterValues(result, comm, tag);
    }

    return fop(src, dest);
}


// ************************************************************************* //
