/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "ccmBase.H"
#include "ccmInternal.H"  // include last to avoid any strange interactions


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::ccm::base::assertNoError
(
    int err,
    const char* msg
)
{
    return ccmGlobalState::assertNoError(static_cast<CCMIOError>(err), msg);
}


bool Foam::ccm::base::assertNoError
(
    int err,
    const std::string& msg
)
{
    return ccmGlobalState::assertNoError(static_cast<CCMIOError>(err), msg);
}


bool Foam::ccm::base::assertNoError(const char* msg) const
{
    return globalState_->assertNoError(msg);
}


bool Foam::ccm::base::assertNoError(const std::string& msg) const
{
    return globalState_->assertNoError(msg);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ccm::base::base()
:
    globalState_(new ccmGlobalState)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ccm::base::~base()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ccm::base::close()
{
    if (globalState_)
    {
        if (CCMIOIsValidEntity(globalState_->root))
        {
            CCMIOCloseFile(nullptr, globalState_->root);
        }
        globalState_.reset(nullptr);

        return true;
    }

    return false;
}


// ************************************************************************* //
