/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::gradingDescriptor

Description
    Handles the specification for grading within a section of a block

    The grading specification is handled is controlled by three parameters:

    - blockFraction: the fraction of the block the section occupies

    - nDivFraction: the fraction of the divisions of the block allocated to
        the section

    - expansionRatio:
       the expansion ratio for the grading with the section of
       block defined as the ratio of end-size / start-size for the section.
       A negative value is trapped and treated as its inverse.

SourceFiles
    gradingDescriptor.C

\*---------------------------------------------------------------------------*/

#ifndef gradingDescriptor_H
#define gradingDescriptor_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class gradingDescriptor;
class gradingDescriptors;

Istream& operator>>(Istream&, gradingDescriptor&);
Ostream& operator<<(Ostream&, const gradingDescriptor&);

/*---------------------------------------------------------------------------*\
                     Class gradingDescriptor Declaration
\*---------------------------------------------------------------------------*/

class gradingDescriptor
{
    // Private Data

        scalar blockFraction_;
        scalar nDivFraction_;
        scalar expansionRatio_;


public:

    // Friendship with gradingDescriptors
    friend class gradingDescriptors;


    // Constructors

        //- Default construct (1, 1, 1)
        gradingDescriptor();

        //- Construct from components
        gradingDescriptor
        (
            const scalar blockFraction,
            const scalar nDivFraction,
            const scalar expansionRatio
        );

        //- Construct from expansionRatio
        explicit gradingDescriptor(const scalar expansionRatio);

        //- Construct from Istream
        explicit gradingDescriptor(Istream& is);


    //- Destructor
    ~gradingDescriptor() = default;


    // Member Functions

        scalar blockFraction() const
        {
            return blockFraction_;
        }

        scalar nDivFraction() const
        {
            return nDivFraction_;
        }

        scalar expansionRatio() const
        {
            return expansionRatio_;
        }

        //- Adjust expansion ratio.
        //  Trap negative value and treat as its inverse.
        void correct();

        //- Return the inverse gradingDescriptor with 1/expansionRatio
        gradingDescriptor inv() const;


    // Member Operators

        bool operator==(const gradingDescriptor&) const;
        bool operator!=(const gradingDescriptor&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, gradingDescriptor&);
        friend Ostream& operator<<(Ostream&, const gradingDescriptor&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
