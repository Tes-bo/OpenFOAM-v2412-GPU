/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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
    Foam::pointHistory

Description

SourceFiles
    pointHistory.C

\*---------------------------------------------------------------------------*/

#ifndef pointHistory_H
#define pointHistory_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointHistory Declaration
\*---------------------------------------------------------------------------*/

class pointHistory
:
    public functionObject
{
    // Private Data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- History point ID
        label historyPointID_;

        //- History point
        vector refHistoryPoint_;

        //- Processor of history point
        label processor_;

        //- Output file name
        word fileName_;

        //- Output file stream
        autoPtr<OFstream> historyFilePtr_;

    // Private Member Functions

        //- Write data
        bool writeData();

        //- No copy construct
        pointHistory(const pointHistory&) = delete;

        //- No copy assignment
        void operator=(const pointHistory&) = delete;

public:

    //- Runtime type information
    TypeName("pointHistory");


    // Constructors

        //- Construct from components
        pointHistory
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    // Member Functions

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- No-op
        virtual bool write()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
