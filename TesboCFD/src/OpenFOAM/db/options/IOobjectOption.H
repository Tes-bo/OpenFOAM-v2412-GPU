/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::IOobjectOption

Description
    A simple container of IOobject preferences.
    Can also be used for general handling of read/no-read/read-if-present
    logic outside of an IOobject.

See also
    Foam::IOobject

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOobjectOption_H
#define Foam_IOobjectOption_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IOobjectOption Declaration
\*---------------------------------------------------------------------------*/

class IOobjectOption
{
public:

    // Public Data Types

        //- Enumeration defining read preferences
        //  Lowest bit encodes 'must read'.
        //  Possible (future) named variants (none | normal | modified | lazy)
        enum readOption : unsigned char
        {
            //! Nothing to be read
            NO_READ = 0,

            //! Reading required
            MUST_READ = 0x1,

            //! Reading required, file watched for runTime modification
            //! [identical to MUST_READ_IF_MODIFIED]
            READ_MODIFIED = 0x3,

            //! Reading required, file watched for runTime modification
            //! [identical to READ_MODIFIED]
            MUST_READ_IF_MODIFIED = 0x3,

            //! Reading is optional [identical to READ_IF_PRESENT]
            LAZY_READ = 0x4,

            //! Reading is optional [identical to LAZY_READ]
            READ_IF_PRESENT = 0x4
        };

        //- Enumeration defining write preferences
        enum writeOption : unsigned char
        {
            //! Ignore writing from objectRegistry::writeObject()
            NO_WRITE = 0,

            //! Automatically write from objectRegistry::writeObject()
            AUTO_WRITE = 0x10
        };

        //- Enumeration for use with registerObject().
        //- Values map to bool (false/true)
        enum registerOption : unsigned char
        {
            //! Do not request registration (bool: false)
            NO_REGISTER = 0,

            //! Request registration (bool: true)
            REGISTER = 1,

            //! Legacy/default registration request (bool: true)
            LEGACY_REGISTER = 2
        };


private:

    // Private Data

        //- Read option
        readOption readOpt_;

        //- Write option
        writeOption writeOpt_;

        //- Should created objects be registered?
        bool registerObject_;

        //- Is object same for all processors?
        bool globalObject_;


public:

    // Constructors

        //- Default construct (NO_READ, NO_WRITE, REGISTER, non-global)
        //- or construct with specified options
        constexpr IOobjectOption
        (
            readOption rOpt = readOption::NO_READ,
            writeOption wOpt = writeOption::NO_WRITE,
            registerOption registerObject = registerOption::REGISTER,
            bool globalObject = false
        ) noexcept
        :
            readOpt_(rOpt),
            writeOpt_(wOpt),
            registerObject_(registerObject),
            globalObject_(globalObject)
        {}

        //- Construct NO_WRITE with specified read/register options
        constexpr IOobjectOption
        (
            readOption rOpt,
            registerOption registerObject = registerOption::REGISTER,
            bool globalObject = false
        ) noexcept
        :
            readOpt_(rOpt),
            writeOpt_(writeOption::NO_WRITE),
            registerObject_(registerObject),
            globalObject_(globalObject)
        {}

        //- Construct NO_READ with specified write/register options
        constexpr IOobjectOption
        (
            writeOption wOpt,
            registerOption registerObject = registerOption::REGISTER,
            bool globalObject = false
        ) noexcept
        :
            readOpt_(readOption::NO_READ),
            writeOpt_(wOpt),
            registerObject_(registerObject),
            globalObject_(globalObject)
        {}

        //- Construct (NO_READ, NO_WRITE) with specified register option
        constexpr IOobjectOption
        (
            registerOption registerObject,
            bool globalObject = false
        ) noexcept
        :
            readOpt_(readOption::NO_READ),
            writeOpt_(writeOption::NO_WRITE),
            registerObject_(registerObject),
            globalObject_(globalObject)
        {}

        //- Construct from components
        //- with specified register option as bool
        constexpr IOobjectOption
        (
            readOption rOpt,
            writeOption wOpt,
            bool registerObject,
            bool globalObject = false
        ) noexcept
        :
            readOpt_(rOpt),
            writeOpt_(wOpt),
            registerObject_(registerObject ? REGISTER : NO_REGISTER),
            globalObject_(globalObject)
        {}

        //- Construct (NO_READ, NO_WRITE)
        //- with specified register option as bool
        explicit constexpr IOobjectOption
        (
            bool registerObject,
            bool globalObject = false
        ) noexcept
        :
            readOpt_(readOption::NO_READ),
            writeOpt_(writeOption::NO_WRITE),
            registerObject_(registerObject ? REGISTER : NO_REGISTER),
            globalObject_(globalObject)
        {}


    // Member Functions

        //- Get the read option
        readOption readOpt() const noexcept { return readOpt_; }

        //- Set the read option \return the previous value
        readOption readOpt(readOption opt) noexcept
        {
            readOption old(readOpt_);
            readOpt_ = opt;
            return old;
        }

        //- Get the write option
        writeOption writeOpt() const noexcept { return writeOpt_; }

        //- Set the write option \return the previous value
        writeOption writeOpt(writeOption opt) noexcept
        {
            writeOption old(writeOpt_);
            writeOpt_ = opt;
            return old;
        }

        //- Should objects created with this IOobject be registered?
        bool registerObject() const noexcept { return registerObject_; }

        //- Change registration preference \return previous value
        bool registerObject(bool on) noexcept
        {
            bool old(registerObject_);
            registerObject_ = on;
            return old;
        }

        //- True if object is treated the same for all processors
        bool globalObject() const noexcept { return globalObject_; }

        //- Change global-object status \return previous value
        bool globalObject(bool on) noexcept
        {
            bool old(globalObject_);
            globalObject_ = on;
            return old;
        }


    // Checks

        //- True if any reading may be required (ie, != NO_READ)
        static bool isAnyRead(readOption opt) noexcept
        {
            return (opt != readOption::NO_READ);
        }

        //- True if any reading may be required (ie, != NO_READ)
        bool isAnyRead() const noexcept
        {
            return (readOpt_ != readOption::NO_READ);
        }

        //- True if (MUST_READ | READ_MODIFIED) bits are set
        static bool isReadRequired(readOption opt) noexcept
        {
            return (opt & readOption::MUST_READ);
        }

        //- True if (MUST_READ | READ_MODIFIED) bits are set
        bool isReadRequired() const noexcept
        {
            return (readOpt_ & readOption::MUST_READ);
        }

        //- True if (LAZY_READ) bits are set [same as READ_IF_PRESENT]
        static bool isReadOptional(readOption opt) noexcept
        {
            return (opt == readOption::LAZY_READ);
        }

        //- True if (LAZY_READ) bits are set [same as READ_IF_PRESENT]
        bool isReadOptional() const noexcept
        {
            return (readOpt_ == readOption::LAZY_READ);
        }

        //- Downgrade readOption optional (LAZY_READ), leaves NO_READ intact.
        static readOption lazierRead(readOption opt) noexcept
        {
            return (opt == readOption::NO_READ ? opt : readOption::LAZY_READ);
        }


    // Housekeeping

        //- Access to the read option
        //  \deprecated(2021-03) - use readOpt(readOption)
        readOption& readOpt() noexcept { return readOpt_; }

        //- Access to the write option
        //  \deprecated(2021-03) - use writeOpt(writeOption)
        writeOption& writeOpt() noexcept { return writeOpt_; }

        //- Access to the register object option
        //  \deprecated(2021-03) - use registerObject(bool)
        bool& registerObject() noexcept { return registerObject_; }

        //- Access to the global object option
        //  \deprecated(2021-03) - use globalObject(bool)
        bool& globalObject() noexcept { return globalObject_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
