/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::InfoProxy

Description
    A helper class for outputting values to Ostream

\*---------------------------------------------------------------------------*/

#ifndef Foam_InfoProxy_H
#define Foam_InfoProxy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Ostream;

template<class T> class InfoProxy;
template<class T> Ostream& operator<<(Ostream&, const InfoProxy<T>&);


/*---------------------------------------------------------------------------*\
                           Class InfoProxy Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class InfoProxy
{
public:

    // Data

        //- The captured object reference
        const T& t_;


    // Constructors

        //- Implicit construct from object reference
        InfoProxy(const T& item) noexcept
        :
            t_(item)
        {}


    // Member Functions

        //- Return const reference to the object
        const T& operator*() const noexcept { return t_; }

        // // TBD: Dereference const pointer to the object
        // const T* operator->() const noexcept { return &t_; }


    // Output Functions

        //- Output information to stream
        friend Ostream& operator<< <T> (Ostream&, const InfoProxy<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
