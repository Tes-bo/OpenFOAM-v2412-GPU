/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::IOMap

Description
    A Map of objects of type \<T\> with automated input and output.
    Is a global object; i.e. can be read from undecomposed case.

SourceFiles
    IOMap.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOMap_H
#define Foam_IOMap_H

#include "Map.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class IOMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOMap
:
    public regIOobject,
    public Map<T>
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool readIOcontents();

public:

    //- The underlying content type
    typedef Map<T> content_type;

    //- Runtime type information
    TypeName("Map");


    // Constructors

        //- Default copy construct
        IOMap(const IOMap&) = default;

        //- Construct from IOobject
        explicit IOMap(const IOobject& io);

        //- Construct from IOobject and size of Map
        IOMap(const IOobject& io, const label size);

        //- Construct from IOobject and a copy of Map content
        IOMap(const IOobject&, const Map<T>& content);

        //- Construct by transferring the Map content
        IOMap(const IOobject&, Map<T>&& content);


    // Factory Methods

        //- Read and return contents. The IOobject will not be registered
        static Map<T> readContents(const IOobject& io);


    //- Destructor
    virtual ~IOMap() = default;


    // Member Functions

        //- This object is global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- The writeData method for regIOobject write operation
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const IOMap<T>& rhs);

        //- Copy or move assignment of entries
        using Map<T>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for IOMap
template<class T>
struct is_globalIOobject<IOMap<T>> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOMap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
