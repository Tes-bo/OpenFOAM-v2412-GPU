/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "Matrix.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"
#include "contiguous.H"
#include "ListPolicy.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(Istream& is)
:
    mRows_(0),
    nCols_(0),
    v_(nullptr)
{
    this->readMatrix(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
bool Foam::Matrix<Form, Type>::readMatrix(Istream& is)
{
    // Anull matrix
    clear();

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck("readMatrix : reading first token");

    if (firstToken.isLabel())
    {
        mRows_ = firstToken.labelToken();
        nCols_ = readLabel(is);
        doAlloc();

        // The total size
        const label len = size();

        if (is.format() == IOstreamOption::BINARY && is_contiguous<Type>::value)
        {
            // Binary and contiguous

            if (len)
            {
                Detail::readContiguous<Type>
                (
                    is,
                    this->data_bytes(),
                    this->size_bytes()
                );

                is.fatalCheck("readMatrix : reading the binary block");
            }
        }
        else
        {
            // Begin of contents marker
            char listDelimiter = is.readBeginList("Matrix");

            if (len)
            {
                if (listDelimiter == token::BEGIN_LIST)
                {
                    auto iter = this->begin();

                    // Loop over rows
                    for (label i = 0; i < mRows_; ++i)
                    {
                        listDelimiter = is.readBeginList("MatrixRow");

                        for (label j = 0; j < nCols_; ++j, (void)++iter)
                        {
                            is >> *iter;
                            is.fatalCheck("readMatrix : reading entry");
                        }

                        is.readEndList("MatrixRow");
                    }
                }
                else  // BEGIN_BLOCK
                {
                    Type elem;
                    is >> elem;

                    is.fatalCheck("readMatrix : reading the single entry");

                    std::fill_n(begin(), size(), elem);
                }
            }

            // End of contents marker
            is.readEndList("Matrix");
        }

        return len;
    }

    FatalIOErrorInFunction(is)
        << "incorrect first token, expected <int>, found "
        << firstToken.info() << nl
        << exit(FatalIOError);

    return 0;
}


template<class Form, class Type>
Foam::Ostream& Foam::Matrix<Form, Type>::writeMatrix
(
    Ostream& os,
    const label shortLen
) const
{
    const Matrix<Form, Type>& mat = *this;
    const label len = mat.size();  // Total size (rows * cols)

    auto iter = mat.cbegin();  // element-wise iterator

    // Rows, columns size
    os  << mat.nRows() << token::SPACE << mat.nCols();

    if (os.format() == IOstreamOption::BINARY && is_contiguous<Type>::value)
    {
        // Binary and contiguous

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write(mat.cdata_bytes(), mat.size_bytes());
        }
    }
    else if (is_contiguous<Type>::value && len > 1 && mat.uniform())
    {
        // Two or more entries, and all entries have identical values.
        os  << token::BEGIN_BLOCK << *iter << token::END_BLOCK;
    }
    else if
    (
        (len <= 1 || !shortLen)
     || (len <= shortLen && is_contiguous<Type>::value)
    )
    {
        // Single-line output (entire matrix)

        // Begin matrix
        os  << token::BEGIN_LIST;

        // Loop over rows
        for (label i = 0; i < mat.nRows(); ++i)
        {
            // Begin row
            os  << token::BEGIN_LIST;

            // Write row
            for (label j = 0; j < mat.nCols(); ++j, (void)++iter)
            {
                if (j) os << token::SPACE;
                os << *iter;
            }

            // End row
            os << token::END_LIST;
        }

        // End matrix
        os << token::END_LIST;
    }
    else if
    (
        (mat.nCols() <= 1 || !shortLen)
     || (mat.nCols() <= shortLen && is_contiguous<Type>::value)
    )
    {
        // Multi-line matrix, single-line rows

        // Begin matrix
        os  << nl << token::BEGIN_LIST;

        // Loop over rows
        for (label i = 0; i < mat.nRows(); ++i)
        {
            // Begin row
            os  << nl << token::BEGIN_LIST;

            // Write row
            for (label j = 0; j < mat.nCols(); ++j, (void)++iter)
            {
                if (j) os << token::SPACE;
                os << *iter;
            }

            // End row
            os << token::END_LIST;
        }

        // End matrix
        os << nl << token::END_LIST << nl;
    }
    else
    {
        // Multi-line output

        // Begin matrix
        os  << nl << token::BEGIN_LIST;

        // Loop over rows
        for (label i=0; i < mat.nRows(); ++i)
        {
            // Begin row
            os  << nl << token::BEGIN_LIST;

            // Write row
            for (label j = 0; j < mat.nCols(); ++j, (void)++iter)
            {
                os << nl << *iter;
            }

            // End row
            os << nl << token::END_LIST;
        }

        // End matrix
        os << nl << token::END_LIST << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
