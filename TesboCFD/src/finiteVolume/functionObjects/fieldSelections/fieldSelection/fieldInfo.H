/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::functionObjects::fieldInfo

Description
    Helper class to store a wordRe and label used by
    Foam::functionObjects::fieldSelection

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_fieldInfo_H
#define Foam_functionObjects_fieldInfo_H

#include "label.H"
#include "wordRes.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward Declarations
class fieldInfo;
Istream& operator>>(Istream&, fieldInfo&);
Ostream& operator<<(Ostream&, const fieldInfo&);

/*---------------------------------------------------------------------------*\
                          Class fieldInfo Declaration
\*---------------------------------------------------------------------------*/

class fieldInfo
{
    // Private Data

        //- Pattern for the field name(s)
        wordRe name_;

        //- Field component
        label component_;

        //- Found the field
        mutable bool found_;


public:

    // Constructors

        //- Default construct
        fieldInfo()
        :
            name_(),
            component_(-1),
            found_(false)
        {}

        //- Construct from components
        explicit fieldInfo(const wordRe& name, const label component = -1)
        :
            name_(name),
            component_(component),
            found_(false)
        {}

        //- Construct from stream
        explicit fieldInfo(Istream& is)
        :
            name_(is),
            component_(readLabel(is)),
            found_(false)
        {}


    //- Destructor
    ~fieldInfo() = default;


    // Member Functions

        //- Return the selector pattern for the field name(s)
        const wordRe& name() const noexcept { return name_; }

        //- Return the component
        label component() const noexcept { return component_; }

        //- Return the found state
        bool found() const noexcept { return found_; }

        //- Set the found state to be 'on'
        void found(bool on) const noexcept { found_ = on; }


        friend bool operator==(const fieldInfo& a, const fieldInfo& b)
        {
            return
            (
                a.found() == b.found()
             && a.component() == b.component()
             && a.name() == b.name()
            );
        }

        friend bool operator!=(const fieldInfo& a, const fieldInfo& b)
        {
            return !(a == b);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, fieldInfo& fi)
        {
            is  >> fi.name_ >> fi.component_ >> fi.found_;
            return is;
        }
        friend Ostream& operator<<(Ostream& os, const fieldInfo& fi)
        {
            os  << fi.name_ << ' '
                << fi.component_ << ' '
                << Switch::name(fi.found_);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
