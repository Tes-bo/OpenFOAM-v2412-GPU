/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::IOPosition

Description
    Helper IO class to read and write particle coordinates (positions).

SourceFiles
    IOPosition.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOPosition_H
#define Foam_IOPosition_H

#include "cloud.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IOPosition Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IOPosition
:
    public regIOobject
{
    // Private Data

        cloud::geometryType geometryType_;

        //- Reference to the cloud
        const CloudType& cloud_;


public:

    // Static Data

        //- Runtime type name information. Use cloud type.
        virtual const word& type() const
        {
            return Cloud<typename CloudType::particleType>::typeName;
        }


    // Constructors

        //- Construct from cloud
        explicit IOPosition
        (
            const CloudType& c,
            cloud::geometryType geomType = cloud::geometryType::COORDINATES
        );


    // Member Functions

        //- Inherit readData from regIOobject
        using regIOobject::readData;

        virtual void readData(Istream&, CloudType&);

        virtual bool write(const bool writeOnProc = true) const;

        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOPosition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
