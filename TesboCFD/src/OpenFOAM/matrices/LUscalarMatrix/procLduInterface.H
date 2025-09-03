/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::procLduInterface

Description
    IO interface for processorLduInterface

SourceFiles
    procLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_procLduInterface_H
#define Foam_procLduInterface_H

#include "labelList.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class lduInterfaceField;
class procLduInterface;

Istream& operator>>(Istream&, procLduInterface&);
Ostream& operator<<(Ostream&, const procLduInterface&);

/*---------------------------------------------------------------------------*\
                      Class procLduInterface Declaration
\*---------------------------------------------------------------------------*/

class procLduInterface
{
    // Private Data

        labelList faceCells_;
        scalarField coeffs_;
        label myProcNo_;
        label neighbProcNo_;
        label tag_;
        label comm_;

public:

    //- Friendship
    friend class LUscalarMatrix;


    // Generated Methods

        //- No copy construct
        procLduInterface(const procLduInterface&) = delete;


    // Constructors

        procLduInterface
        (
            const lduInterfaceField& interface,
            const scalarField& coeffs
        );

        //- Read construct from Istream
        explicit procLduInterface(Istream& is);

        autoPtr<procLduInterface> clone()
        {
            NotImplemented;
            return nullptr;
        }

        static autoPtr<procLduInterface> New(Istream& is)
        {
            return autoPtr<procLduInterface>::New(is);
        }


    // IO Operations

        void read(Istream& is);
        void write(Ostream& os) const;

        friend Istream& operator>>(Istream&, procLduInterface&);
        friend Ostream& operator<<(Ostream&, const procLduInterface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
