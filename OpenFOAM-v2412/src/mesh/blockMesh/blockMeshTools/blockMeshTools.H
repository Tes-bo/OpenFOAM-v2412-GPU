/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

Namespace
    Foam::blockMeshTools

Description
    Tools for parsing label(List) with dictionary lookup.

SourceFiles
    blockMeshTools.C
    blockMeshToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_blockMeshTools_H
#define Foam_blockMeshTools_H

#include "dictionary.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockMeshTools
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- In-place read with dictionary lookup
    void read(Istream&, label& val, const dictionary&);

    //- In-place read with dictionary lookup
    template<class T>
    void read(Istream&, List<T>& list, const dictionary&);

    //- Return-read with dictionary lookup
    label read(Istream&, const dictionary&);

    //- Return-read with dictionary lookup
    template<class T>
    List<T> read(Istream& is, const dictionary&);

    //- Write with dictionary lookup
    void write(Ostream&, const label val, const dictionary&);

    //- Linear search for labelled entry, nullptr if not found
    const entry* findEntry(const dictionary& dict, const label val);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockMeshTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blockMeshToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
