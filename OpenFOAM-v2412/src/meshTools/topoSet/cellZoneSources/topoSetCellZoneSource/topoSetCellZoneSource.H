/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::topoSetCellZoneSource

Description
    The \c topoSetCellZoneSource is a intermediate class
    for handling \c topoSet sources for selecting cell zones.

    The class is not an executable \c topoSet itself,
    yet a provider for common entries to its derived \c topoSet.

SourceFiles
    topoSetCellZoneSource.C

\*---------------------------------------------------------------------------*/

#ifndef topoSetCellZoneSource_H
#define topoSetCellZoneSource_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellZoneSet;

/*---------------------------------------------------------------------------*\
                    Class topoSetCellZoneSource Declaration
\*---------------------------------------------------------------------------*/

class topoSetCellZoneSource
:
    public topoSetSource
{
public:

    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetCellZoneSource,
            word,
            (
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

        // For the Istream constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetCellZoneSource,
            istream,
            (
                const polyMesh& mesh,
                Istream& is
            ),
            (mesh, is)
        );


    // Constructors

        //- Construct from mesh
        explicit topoSetCellZoneSource(const polyMesh& mesh);

        //- Construct from mesh, use "verbose" entry if present
        topoSetCellZoneSource(const polyMesh& mesh, const dictionary& dict);


        //- Clone (disallowed)
        autoPtr<topoSetCellZoneSource> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a reference to the selected source type
        static autoPtr<topoSetCellZoneSource> New
        (
            const word& sourceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Return a reference to the selected source type
        static autoPtr<topoSetCellZoneSource> New
        (
            const word& sourceType,
            const polyMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~topoSetCellZoneSource() = default;


    // Member Functions

        //- The source category is a cellZone
        virtual topoSetSource::sourceType setType() const
        {
            return CELLZONE_SOURCE;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
