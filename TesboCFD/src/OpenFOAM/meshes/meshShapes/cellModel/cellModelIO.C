/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "cellModel.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellModel::cellModel(Istream& is)
{
    const dictionaryEntry dictEntry(dictionary::null, is);
    const dictionary& dict = dictEntry.dict();

    name_ = dictEntry.keyword();
    dict.readEntry("index", index_);
    dict.readEntry("numberOfPoints", nPoints_);
    dict.readEntry("faces", faces_);
    dict.readEntry("edges", edges_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const cellModel& cm)
{
    os  << "name" << tab << cm.name() << tab
        << "index" << tab << cm.index() << tab
        << "numberOfPoints" << tab << cm.nPoints() << tab
        << "faces" << tab << cm.modelFaces() << tab
        << "edges" << tab << cm.modelEdges() << endl;

    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<cellModel>& iproxy
)
{
    const auto& cm = *iproxy;

    os  << "name = " << cm.name() << ", "
        << "index = " << cm.index() << ", "
        << "number of points = " << cm.nPoints() << ", "
        << "number of faces = " << cm.nFaces() << ", "
        << "number of edges = " << cm.nEdges() << endl;

    return os;
}


// ************************************************************************* //
