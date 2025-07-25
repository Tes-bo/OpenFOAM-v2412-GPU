/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::Function1Types::TableFile

Description
    Templated table container function where data are read from file.

    Usage:
    \verbatim
        <entryName> tableFile;
        <entryName>Coeffs
        {
            file                dataFile;    // name of data file
            outOfBounds         clamp;       // optional out-of-bounds handling
            interpolationScheme linear;      // optional interpolation method
        }
    \endverbatim

    Data are stored as a list of Tuple2's. First column is always stored as
    scalar entries.  Data are read in the form, e.g. for an entry \<entryName\>
    that is (scalar, vector):
    \verbatim
        (
            (0.0 (1 2 3))
            (1.0 (4 5 6))
        );
    \endverbatim

See Also
    Foam::Function1Types::Table

SourceFiles
    TableFile.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_TableFile_H
#define Foam_Function1Types_TableFile_H

#include "TableBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class TableFile Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableFile
:
    public TableBase<Type>
{
    // Private Data

        //- File name for table (mandatory)
        fileName fName_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const TableFile<Type>&) = delete;


public:

    //- Runtime type information
    TypeName("tableFile");


    // Constructors

        //- Construct from entry name and "file" found in dictionary
        TableFile
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        explicit TableFile(const TableFile<Type>& tbl);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~TableFile() = default;


    // Member Functions

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableFile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
