/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    Foam::Gather

Description
    Gather data from all processors onto all processors.

SourceFiles
    Gather.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_TestGather_H
#define Foam_TestGather_H

#include "List.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Gather Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Gather
:
    public List<Type>
{
    // Private Data

        //- Number of processors (1 for sequential)
        label nProcs_;

        //- Storage
        //List<Type> list_;


public:

    // Constructors

        //- Construct from containers across processors
        Gather(const Type& localData, const bool redistribute=true);


    // Member Functions

        List<Type>& list() noexcept { return *this; }
        const List<Type>& list() const noexcept { return *this; }
    //  List<Type>& list() noexcept { return list_; }
    //  const List<Type>& list() const noexcept { return list_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Gather.txx"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
