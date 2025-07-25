/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
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

Class
    Foam::badQualityToCell

Description
    Selects bad quality cells (using snappyHexMesh/cvMesh mesh quality selector)

SourceFiles
    badQualityToCell.C

\*---------------------------------------------------------------------------*/

#ifndef badQualityToCell_H
#define badQualityToCell_H

#include "topoSetCellSource.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class badQualityToCell Declaration
\*---------------------------------------------------------------------------*/

class badQualityToCell
:
    public topoSetCellSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Mesh quality dictionary
        const dictionary dict_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("badQualityToCell");

    // Constructors

        //- Construct from dictionary
        badQualityToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        badQualityToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~badQualityToCell() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
