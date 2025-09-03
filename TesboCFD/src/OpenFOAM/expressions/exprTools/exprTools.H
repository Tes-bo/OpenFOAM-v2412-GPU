/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Namespace
    Foam::exprTools

Description
    Helper methods for creating expressions.

SourceFiles
    exprTools.C
    exprToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_exprTools_H
#define Foam_exprTools_H

#include "exprOps.H"
#include "exprString.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Namespace exprTools Declarations
\*---------------------------------------------------------------------------*/

namespace exprTools
{

//- Stringified value for zero of given type
template<class Type>
string zeroValue();

//- Stringified version of data with comma separators between components
template<class Type>
string toString(const Type& data, const word& prefix);

//- Stringified version of data with comma separators between components.
//- Uses prefix corresponding to the pTraits of Type.
template<class Type>
string toString(const Type& data);

//- Stringified version of primitiveEntry of Type.
//- Uses comma separators and prefix corresponding to the pTraits of Type
template<class Type>
string toString(ITstream& is);

//- Get an expression string list from a dictionary
//  Supports list and inline specifications, as well as indirect
//  inclusion of other expression string lists.
//
//  Syntax examples,
//  \code
//
//      // Single items. Inline or list format.
//      list1a  "var1=p";
//      list1b  ( "var1=p" );
//
//      // Multiple items. Inline or list format.
//      list2a  "var2=T; var3=U; var4=rho";
//      list2b  ( "var2=T" "var3=U; var4=rho" );
//
//      // Indirection. Inline or list format.
//      list3a  "#list1a; #list2b";
//      list2b  ( "var1=p" "#list2a" );
//
//  \endcode
List<expressions::exprString> getList
(
    const word& keyword,
    const dictionary& dict,
    bool mandatory = true,
    label recursionDepth = 0
);


} // End namespace exprTools


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "exprToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
