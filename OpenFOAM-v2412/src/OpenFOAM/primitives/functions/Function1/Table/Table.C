/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "Table.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Table<Type>::Table
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    TableBase<Type>(entryName, dict, obrPtr),
    fName_()
{
    const entry* eptr = dict.findEntry(entryName, keyType::LITERAL);

    if (eptr && eptr->isStream())
    {
        // Primitive (inline) format. Eg,
        // key table ((0 0) (10 1));

        ITstream& is = eptr->stream();
        if (is.peek().isWord())
        {
            is.skip();  // Discard leading 'table'
        }
        is >> this->table_;
        dict.checkITstream(is, entryName);
    }
    else if (dict.readIfPresent("file", fName_))
    {
        // Dictionary format - "file" lookup. Eg,
        // key { type table; file "name"; }

        fileName expandedFile(fName_);
        expandedFile.expand();

        autoPtr<ISstream> isPtr(fileHandler().NewIFstream(expandedFile));
        if (isPtr && isPtr->good())
        {
            *isPtr >> this->table_;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Cannot open file: " << expandedFile << nl
                << exit(FatalIOError);
        }
    }
    else
    {
        // Dictionary format - "values" lookup. Eg,
        //
        // key { type table; values ((0 0) (10 1)); }

        dict.readEntry("values", this->table_);
    }

    TableBase<Type>::initialise();
}


template<class Type>
Foam::Function1Types::Table<Type>::Table(const Table<Type>& tbl)
:
    TableBase<Type>(tbl),
    fName_(tbl.fName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Table<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));

    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    if (fName_.empty())
    {
        os.writeEntry("values", this->table_);
    }
    else
    {
        os.writeEntry("file", fName_);
    }

    os.endBlock();
}


// ************************************************************************* //
