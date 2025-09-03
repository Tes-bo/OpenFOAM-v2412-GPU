/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::fileFormats::surfaceFormatsCore

Description
    A collection of helper functions for reading/writing surface formats.

SourceFiles
    surfaceFormatsCore.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFormatsCore_H
#define surfaceFormatsCore_H

#include "Map.H"
#include "HashSet.H"
#include "wordRes.H"
#include "labelList.H"
#include "surfZoneList.H"
#include "surfZoneIdentifierList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ISstream;
class Time;
class triSurface;
class IOobject;
class dictionary;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
               Class fileFormats::surfaceFormatsCore Declaration
\*---------------------------------------------------------------------------*/

class surfaceFormatsCore
{
protected:

    // Friendship with triSurface (for IO helpers)
    friend class ::Foam::triSurface;


    // General

        //- Read non-empty and non-comment line
        static string getLineNoComment(ISstream& is, const char comment='#');

        //- Return a surfZone list with a single entry, the size of which
        //- corresponds to that of the container
        template<class Container>
        static List<surfZone> oneZone
        (
            const Container& container,
            const word& name = "zone0"
        )
        {
            return List<surfZone>(1, surfZone(name, container.size()));
        }

        //- Return ids for zone/patch that match by name.
        //  Uses a combination of allow and deny lists.
        //
        //  \see Foam::stringListOps::findMatching for details about matching
        static labelList getSelectedPatches
        (
            const surfZoneList& patches,
            const wordRes& allow,
            const wordRes& deny = wordRes()
        );


    // IO helpers

        //- Return fileName.
        //  If fileName is relative gets treated as local to IOobject.
        //  Duplicate of triSurfaceMesh::relativeFilePath
        static fileName relativeFilePath
        (
            const IOobject& io,
            const fileName& f,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Return fileName to load IOobject from.
        //  Fatal if the file does not exist
        static fileName checkFile
        (
            const IOobject& io,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Return fileName to load IOobject from.
        //  Supports optional override of fileName with "file" entry
        //  Fatal if the file does not exist
        static fileName checkFile
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );


public:

    // Static Data

        //- The file extension corresponding to 'native' surface format
        //  Normally "ofs" (mnemonic: OF = OpenFOAM, S = Surface)
        static word nativeExt;


    // Static Member Functions

        //- Verbose checking of fileType in the list of available types
        static bool checkSupport
        (
            const wordHashSet& available,
            const word& fileType,
            const bool verbose = false,
            const char* functionName = nullptr
        );

        //- Use IOobject information to resolve file to load from,
        //- or empty if the file does not exist.
        static fileName findFile
        (
            const IOobject& io,
            const bool isGlobal = true  //!< resolve as a global file
        );

        //- Use IOobject information to resolve file to load from,
        //- or empty if the file does not exist.
        //  Supports optional override of fileName with "file" entry
        static fileName findFile
        (
            const IOobject& io,
            const dictionary& dict,
            const bool isGlobal = true  //!< resolve as a global file
        );


//        //- Return the local file name (within time directory)
//        //  NEEDS FIXING
//        static fileName localMeshFileName
//        (
//            const word& surfName = word::null
//        );
//
//        //- Find instance with surfName
//        //  NEEDS FIXING
//        static fileName findMeshInstance
//        (
//            const Time&,
//            const word& surfName = word::null
//        );
//
//        //- Find mesh file with surfName
//        //  NEEDS FIXING
//        static fileName findMeshFile
//        (
//            const Time&,
//            const word& surfName = word::null
//        );


    // Constructors

        //- Default construct
        surfaceFormatsCore() = default;


    //- Destructor
    virtual ~surfaceFormatsCore() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
