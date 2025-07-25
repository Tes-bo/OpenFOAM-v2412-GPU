/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::gradingDescriptors

Description
    List of gradingDescriptor for the sections of a block with additional IO
    functionality

SourceFiles
    gradingDescriptors.C

\*---------------------------------------------------------------------------*/

#ifndef gradingDescriptors_H
#define gradingDescriptors_H

#include "gradingDescriptor.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class gradingDescriptors;
Istream& operator>>(Istream& is, gradingDescriptors& gd);

/*---------------------------------------------------------------------------*\
                     Class gradingDescriptors Declaration
\*---------------------------------------------------------------------------*/

class gradingDescriptors
:
    public List<gradingDescriptor>
{
public:

    // Constructors

        //- Default construct with a single default gradingDescriptor
        gradingDescriptors();

        //- Construct with specified number of default gradingDescriptor
        explicit gradingDescriptors(const label len);

        //- Construct from a single gradingDescriptor
        explicit gradingDescriptors(const gradingDescriptor& gd);


    // Member Functions

        //- Adjust expansion ratios.
        //  Trap negative value and treat as its inverse.
        void correct();

        //- Normalize blockFractions and nDivFractions for the list
        //- of gradingDescriptors, and call correct()
        void normalise();

        //- Return the inverse gradingDescriptors with 1/expansionRatio
        gradingDescriptors inv() const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, gradingDescriptors&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
