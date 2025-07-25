/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Simple test program for compiling/running openmp

\*---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <iostream>

#if _OPENMP
#include <omp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::cout
        << "c++ = " << __cplusplus << '\n';

#if _OPENMP
    std::cout
        << "openmp = " << _OPENMP << '\n'
        << "Initial threads = " << omp_get_num_threads() << "\n\n";

    // Fork threads with their own copies of variables
    int nThreads, threadId;

    #pragma omp parallel private(nThreads, threadId)
    {
        threadId = omp_get_thread_num();
        nThreads = omp_get_num_threads();

        // Master thread
        if (threadId == 0)
        {
            // Printf rather than cout to ensure that it emits in one go
            printf("Number of threads = %d\n", nThreads);
        }

        // Printf rather than cout to ensure that it emits in one go
        printf("Called from thread = %d\n", threadId);
    }
#else
    std::cout << "Compiled without openmp!\n";
#endif

    return 0;
}


// ************************************************************************* //
