/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

Description
    Internal bits for wrapping libccmio - do not use directly

\*---------------------------------------------------------------------------*/

#ifndef ccmInternal_H
#define ccmInternal_H

#include "className.H"
#include "List.H"

// Headers and definitions for using libccmio
#include "libccmio/ccmio.h"

// Low-level routines are also needed
#include "libccmio/ccmiocore.h"
#include "libccmio/ccmioutility.h"

// Do not allow these defines to leak through
#undef TRUE
#undef FALSE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// @cond
// internal use only - skip doxygen documentation
namespace Foam
{
namespace ccm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Declare namespace and debug information.
NamespaceName("ccm");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- A C++ wrapper for the C struct
class ccmID : public CCMIOID {};

//- A C++ wrapper for the C struct
class ccmNODE
:
    public CCMIONode
{
public:

    //- Recast constructor
    ccmNODE(CCMIONode node)
    :
       CCMIONode(node)
    {}
};


//- A C++ wrapper for the enum
class ccmDimension
{
    CCMIODimensionality dims_;

public:

    //- Construct from components
    ccmDimension(CCMIODimensionality dims)
    :
        dims_(dims)
    {}

    //- Return underlying enum
    CCMIODimensionality operator()() const
    {
        return dims_;
    }
};


//- Maintain overall global states (error, rootNode)
class ccmGlobalState
{
public:

    //- Maintain error state between calls
    CCMIOError error;

    //- Root node in the CCM file
    ccmID root;

    //- Null constructor. Start with no error, but root in invalid state.
    ccmGlobalState();

    //- True if there is an error
    bool hasError() const
    {
        return (error != kCCMIONoErr);
    }

    //- Die with msg if there is an error
    //  Return true if there is no error
    static bool assertNoError(CCMIOError err, const char *msg);

    //- Die with msg if there is an error
    //  Return true if there is no error
    static bool assertNoError(CCMIOError err, const std::string& msg);

    //- Die with msg if there is an error
    //  Return true if there is no error
    bool assertNoError(const char *msg) const;

    //- Die with msg if there is an error
    //  Return true if there is no error
    bool assertNoError(const std::string& msg) const;

};


//- MapIds for various parts
class ccmMaps
{
public:

    //- The cell map
    ccmID cells;

    //- The internalFaces map
    ccmID internalFaces;

    //- Boundary region maps
    List<ccmID> boundary;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// @endcond
// internal use only - skip doxygen documentation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
