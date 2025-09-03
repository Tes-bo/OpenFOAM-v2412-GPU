/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::logFiles

Description
    functionObject base class for creating, maintaining and writing log
    files e.g. integrated or averaged field data vs time.

See also
    Foam::functionObject
    Foam::functionObjects::writeFile

SourceFiles
    logFiles.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_logFiles_H
#define functionObjects_logFiles_H

#include "writeFile.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class functionObjects::logFiles Declaration
\*---------------------------------------------------------------------------*/

class logFiles
:
    public writeFile
{

protected:

    // Protected data

        //- File names
        wordList names_;

        //- File pointer
        PtrList<OFstream> filePtrs_;


    // Protected Member Functions

        //- Create the output file
        virtual void createFiles();

        //- Reset the list of names from a wordList
        virtual void resetNames(const wordList& names);

        //- Reset the list of names to a single name entry
        virtual void resetName(const word& name);


private:

    // Private Member Functions

        //- No copy construct
        logFiles(const logFiles&) = delete;

        //- No copy assignment
        void operator=(const logFiles&) = delete;


public:

    // Constructors

        //- Construct from objectRegistry and prefix
        logFiles
        (
            const objectRegistry& obr,
            const word& prefix
        );


        //- Construct from objectRegistry and prefix, and read options
        //  from dictionary
        logFiles
        (
            const objectRegistry& obr,
            const word& prefix,
            const dictionary& dict
        );


    //- Destructor
    virtual ~logFiles() = default;


    // Member Functions

        //- Return const access to the names
        const wordList& names() const;

        //- Return access to the files
        PtrList<OFstream>& files();

        //- Return file 'i'
        OFstream& files(const label i);

        //- Write function
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
