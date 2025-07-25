/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "Gather.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Gather<Type>::Gather(const Type& localData, const bool redistribute)
:
    nProcs_(Foam::max(1, UPstream::nProcs()))
{
    this->list().resize(nProcs_);

    //
    // Collect sizes on all processor
    //

    if (UPstream::parRun())
    {
        if (UPstream::master())
        {
            auto iter = this->list().begin();
            *iter = localData;

            // Receive data
            for (const int proci : UPstream::subProcs())
            {
                ++iter;
                IPstream::recv(*iter, proci);
            }

            // Send data
            for (const int proci : UPstream::subProcs())
            {
                if (redistribute)
                {
                    OPstream::send(*this, proci);
                }
                else
                {
                    // Dummy send (to balance sends/receives)
                    OPstream::send(label(0), proci);
                }
            }
        }
        else
        {
            // Send my local data to master
            OPstream::send(localData, UPstream::masterNo());

            // Receive data from master
            if (redistribute)
            {
                IPstream::recv(*this, UPstream::masterNo());
            }
            else
            {
                // Dummy receive
                label dummy;
                IPstream::recv(dummy, UPstream::masterNo());
            }
        }
    }
    else
    {
        this->list().resize(1);
        this->list()[0] = localData;
    }
}


// ************************************************************************* //
