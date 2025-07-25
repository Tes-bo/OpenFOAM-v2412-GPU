/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::binModel::decomposePatchValues
(
    List<List<Type>>& data,
    const label bini,
    const Type& v,
    const vector& n
) const
{
    return decomposePatchValues_;
}


template<class Type>
Foam::string Foam::binModel::writeComponents(const word& stem) const
{
    if (pTraits<Type>::nComponents == 1)
    {
        return stem;
    }

    string result = "";
    for (label cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        if (cmpt) result += " ";
        result += stem + "_" + word(pTraits<Type>::componentNames[cmpt]);
    }
    return result;
};


template<class Type>
void Foam::binModel::writeBinnedData
(
    List<List<Type>>& data,
    Ostream& os
) const
{
    if (cumulative_)
    {
        for (auto& datai : data)
        {
            for (label bini = 1; bini < nBin_; ++bini)
            {
                datai[bini] += datai[bini-1];
            }
        }
    }

    writeCurrentTime(os);

    for (label bini = 0; bini < nBin_; ++bini)
    {
        Type total = Zero;

        for (label i = 0; i < data.size(); ++i)
        {
            total += data[i][bini];
        }

        writeValue(os, total);

        for (label i = 0; i < data.size(); ++i)
        {
            writeValue(os, data[i][bini]);
        }
    }

    os  << endl;
}


// ************************************************************************* //
