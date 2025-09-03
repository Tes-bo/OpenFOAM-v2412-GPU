/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::Function1Types::CSV

Description
    Templated CSV function.

    Reference column is always a scalar, e.g. time.

    Usage:
    \verbatim
        <entryName> csvFile;
        <entryName>Coeffs
        {
            nHeaderLine         4;          // number of header lines
            refColumn           0;          // reference column index
            componentColumns    (1 2 3);    // component column indices
            separator           ",";        // optional (defaults to ",")
            mergeSeparators     no;         // merge multiple separators
            file                "fileXYZ";  // name of csv data file
            outOfBounds         clamp;      // optional out-of-bounds handling
            interpolationScheme linear;     // optional interpolation scheme
        }
    \endverbatim

SourceFiles
    CSV.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_CSV_H
#define Foam_Function1Types_CSV_H

#include "Function1.H"
#include "TableBase.H"
#include "Tuple2.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class CSV Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CSV
:
    public TableBase<Type>
{
    // Private Data

        //- Number header lines
        const label nHeaderLine_;

        //- Column of the time
        const label refColumn_;

        //- Labels of the components
        const labelList componentColumns_;

        //- Separator character
        const char separator_;

        //- Merge separators flag, e.g. ',,,' becomes ','
        const bool mergeSeparators_;

        //- File name for csv table
        fileName fName_;


    // Private Member Functions

        //- Get component columns entry
        static labelList getComponentColumns
        (
            const word& name,
            const dictionary& dict
        );

        //- Read csv data table
        void read();

        //- Read component values from the split string
        Type readValue(const UList<string>& strings) const;


public:

    // Generated Methods

        //- No copy assignment
        void operator=(const CSV<Type>&) = delete;


    //- Declare type-name, virtual type (with debug switch)
    TypeName("csvFile");


    // Constructors

        //- Construct from entry name, dictionary and optional registry
        CSV
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr,
            const fileName& fName = fileName::null
        );

        //- Copy construct
        explicit CSV(const CSV<Type>& csv);

        //- Return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return Function1<Type>::Clone(*this);
        }


    //- Destructor
    virtual ~CSV() = default;


    // Member Functions

        //- Return const access to the file name
        virtual const fileName& fName() const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CSV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
