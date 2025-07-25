/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

Description
    File-local code for setting/resetting signal handlers.

SourceFiles
    signalMacros.C

\*---------------------------------------------------------------------------*/

#include "error.H"
#include <csignal>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Saved old signal trapping setting (file-local variable)
static __p_sig_fn_t oldAction_ = SIG_DFL;


static void resetHandler(const char *what, int sigNum)
{
    const __p_sig_fn_t prev = ::signal(sigNum, oldAction_);
    oldAction_ = SIG_DFL;

    if (SIG_ERR == prev)
    {
        FatalError
            << "Cannot unset " << what << " signal (" << sigNum
            << ") trapping" << endl
            << abort(FatalError);
    }
}


static void setHandler(const char *what, int sigNum, void (*handler)(int))
{
    oldAction_ = ::signal(sigNum, handler);

    if (SIG_ERR == oldAction_)
    {
        FatalError
            << "Could not set " << what << " signal (" << sigNum
            << ") trapping" << endl
            << abort(FatalError);
    }
}

} // End namespace Foam


// ************************************************************************* //
