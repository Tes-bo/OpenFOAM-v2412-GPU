/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::IjkField<Type>::resize
(
    const labelVector& newSizes,
    const Type& val
)
{
    labelVector& ourSizes = sizes();

    if (ijk_.empty() || !cmptProduct(newSizes))
    {
        // Either/both are empty, can redimension directly
        ourSizes = newSizes;
        Field<Type>::resize(ijk_.size(), val);
        return;
    }

    const unsigned diffs
    (
        (ourSizes.x() != newSizes.x() ? 0x100 : 0)
      | (ourSizes.y() != newSizes.y() ? 0x010 : 0)
      | (ourSizes.z() != newSizes.z() ? 0x001 : 0)
    );

    switch (diffs)
    {
        case 0x000:
            // No change
            return;
            break;

        case 0x001:
            // Change in k only, can redimension directly
            ourSizes = newSizes;
            Field<Type>::resize(ijk_.size(), val);
            return;
            break;

        case 0x010:
            // 2D change in j only, can redimension directly
            if (ourSizes.z() == 1)
            {
                ourSizes = newSizes;
                Field<Type>::resize(ijk_.size(), val);
                return;
            }
            break;
    }


    if ((ourSizes.x()*ourSizes.y()) == newSizes.x()*newSizes.y())
    {
        // Re-partition i,j with the same size
        ourSizes = newSizes;
        Field<Type>::resize(ijk_.size(), val);
        return;
    }


    IjkField<Type>& ourContent = *this;

    IjkField<Type> newContent(newSizes, val);

    // Need new addressing space
    const label ni = min(ourSizes.x(), newSizes.x());
    const label nj = min(ourSizes.y(), newSizes.y());
    const label nk = min(ourSizes.z(), newSizes.z());

    for (label k=0; k<nk; ++k)
    {
        for (label j=0; j<nj; ++j)
        {
            for (label i=0; i<ni; ++i)
            {
                newContent(i,j,k) = ourContent(i,j,k);
            }
        }
    }

    ourSizes = newSizes;
    Field<Type>::transfer(newContent);
}


template<class Type>
void Foam::IjkField<Type>::resize(const labelVector& newSizes)
{
    resize(newSizes, Zero);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::IjkField<Type>::operator=(const IjkField<Type>& rhs)
{
    if (this != &rhs)
    {
        sizes() = rhs.sizes();

        Field<Type>::operator=(rhs);
    }
}


template<class Type>
void Foam::IjkField<Type>::operator=(const tmp<IjkField<Type>>& rhs)
{
    if (this != &(rhs()))
    {
        sizes() = rhs().sizes();

        Field<Type>::operator=(rhs());
    }
}


// ************************************************************************* //
