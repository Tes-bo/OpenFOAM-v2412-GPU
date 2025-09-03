/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::solverTemplate

Description
    Class to store solver template specifications

\*---------------------------------------------------------------------------*/

#ifndef solverTemplate_H
#define solverTemplate_H

#include "boolList.H"
#include "wordList.H"
#include "dimensionSet.H"
#include "IOobject.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                       Class solverTemplate Declaration
\*---------------------------------------------------------------------------*/

class solverTemplate
{
public:

    // Public enumerations

        //- Solver type
        enum solverType
        {
            stCompressible,
            stIncompressible,
            stBuoyant,
            stUnknown
        };

        //- Solver type names
        static const Enum<solverType> solverTypeNames_;


private:

    // Private data

        //- Solver type
        solverType solverType_;

        //- Multi-region flag
        bool multiRegion_;


        // Per-region info

            //- Region types
            wordList regionTypes_;

            //- Region names
            wordList regionNames_;

            //- Field names
            List<wordList> fieldNames_;

            //- Field types
            List<wordList> fieldTypes_;

            //- Field dimensions
            List<PtrList<dimensionSet>> fieldDimensions_;


    // Public member functions

        //- Read a word from a dictionary (offers some protection...)
        word readFromDict
        (
            IOobject& dictHeader,
            const word& entryName
        ) const;

        //- Read fluid region templates
        dictionary readFluidFieldTemplates
        (
            const word& regionName,
            const fileName& baseDir,
            const dictionary& solverDict,
            const Time& runTime
        ) const;

        //- Read solid region templates
        dictionary readSolidFieldTemplates
        (
            const word& regionName,
            const dictionary& solverDict
        ) const;

        //- Set the properties for region with index regionI
        void setRegionProperties
        (
            const dictionary& dict,
            const word& regionType,
            const word& regionName,
            const label regionI
        );


public:

    //- Constructor
    solverTemplate
    (
        const fileName& baseDir,
        const Time& runTime,
        const word& regionName
    );


    // Public member functions

        //- Solver type name
        word type() const;

        //- Return the multi-region flag
        bool multiRegion() const;

        //- Return the number of regions
        label nRegion() const;


        // Per-region info

            //- Return the region type
            const word& regionType(const label regionI) const;

            //- Return the region name
            const word& regionName(const label regionI) const;

            //- Return the field names
            const wordList& fieldNames(const label regionI) const;

            //- Return the field types
            const wordList& fieldTypes(const label regionI) const;

            //- Return the field dimensions
            const PtrList<dimensionSet>& fieldDimensions
            (
                const label regionI
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
