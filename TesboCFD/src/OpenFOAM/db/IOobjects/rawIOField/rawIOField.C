/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "rawIOField.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::rawIOField<Type>::readIOcontents
(
    Istream& is,
    IOobjectOption::readOption readAverage
)
{
    is >> static_cast<Field<Type>&>(*this);

    if (IOobjectOption::isReadRequired(readAverage))
    {
        is >> average_;
        hasAverage_ = true;
    }
    else if (IOobjectOption::isReadOptional(readAverage))
    {
        // Slightly heavy handed
        const bool oldThrowingIOerr = FatalIOError.throwing(true);

        try
        {
            is >> average_;
            hasAverage_ = true;
        }
        catch (const Foam::IOerror& err)
        {
            average_ = Zero;
            hasAverage_ = false;
        }
        FatalIOError.throwing(oldThrowingIOerr);
    }
}


template<class Type>
bool Foam::rawIOField<Type>::readIOcontents
(
    IOobjectOption::readOption readAverage
)
{
    if (isReadRequired() || isReadOptional())
    {
        bool haveFile = false;
        bool haveHeader = false;

        // Replacement of regIOobject::headerOk() since that one complains
        // if there is no header. TBD - Move up to headerOk()/fileHandler.
        {
            const fileName fName(filePath());

            // Try to open raw first
            autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fName));

            if (isPtr && isPtr->good())
            {
                haveFile = true;

                auto& is = *isPtr;

                const token firstToken(is);

                haveHeader = is.good() && firstToken.isWord("FoamFile");
            }

            if (debug)
            {
                Pout<< "rawIOField : object:" << name()
                    << " haveFile:" << haveFile
                    << " haveHeader:" << haveHeader << endl;
            }
        }


        if (haveHeader)
        {
            // Read but don't fail upon wrong class. Could extend by providing
            // wanted typeName. Tbd.
            Istream& is = readStream(word::null);

            if (is.good())
            {
                readIOcontents(is, readAverage);
                close();
            }
        }
        else if (haveFile)
        {
            // Failed reading header - fall back to IFstream
            autoPtr<ISstream> isPtr(fileHandler().NewIFstream(objectPath()));

            if (isPtr && isPtr->good())
            {
                readIOcontents(*isPtr, readAverage);
            }
            else
            {
                // Error if required but missing
                if (isReadRequired())
                {
                    FatalIOErrorInFunction(*isPtr)
                        << "Trying to read raw field" << endl
                        << exit(FatalIOError);
                }
            }
        }

        if (debug)
        {
            Pout<< "rawIOField : object:" << name()
                << " size:" << this->size() << endl;
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rawIOField<Type>::rawIOField
(
    const IOobject& io,
    IOobjectOption::readOption readAverage
)
:
    regIOobject(io),
    hasAverage_(false),
    average_(Zero)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<rawIOField<Type>>();

    readIOcontents(readAverage);
}


template<class Type>
Foam::rawIOField<Type>::rawIOField
(
    const IOobject& io,
    const bool readAverage
)
:
    rawIOField<Type>
    (
        io,
        (
            readAverage
          ? IOobjectOption::readOption::MUST_READ
          : IOobjectOption::readOption::NO_READ
        )
    )
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type> Foam::rawIOField<Type>::readContents(const IOobject& io)
{
    IOobject rio(io, IOobjectOption::NO_REGISTER);
    if (rio.readOpt() == IOobjectOption::READ_MODIFIED)
    {
        rio.readOpt(IOobjectOption::MUST_READ);
    }

    rawIOField<Type> reader(rio);

    return Field<Type>(std::move(static_cast<Field<Type>&>(reader)));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::rawIOField<Type>::writeData(Ostream& os) const
{
    os  << static_cast<const Field<Type>&>(*this);
    if (hasAverage_)
    {
        os << token::NL << average_;
    }
    return os.good();
}


// ************************************************************************* //
