/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::DESModelBase

Description
    Base class for DES models providing an interfaces to DES fields.

SourceFiles
    DESModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_DESModelBase_H
#define Foam_DESModelBase_H

#include "className.H"
#include "tmp.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class DESModelBase Declaration
\*---------------------------------------------------------------------------*/

class DESModelBase
{
public:

    //- Runtime type information
    ClassName("DESModelBase");

    // Constructors

        //- Constructor
        DESModelBase() noexcept = default;


    //- Destructor
    virtual ~DESModelBase() = default;


    // Public Member Functions

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const = 0;

        //- Return the shielding function
        virtual tmp<volScalarField> fd() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
