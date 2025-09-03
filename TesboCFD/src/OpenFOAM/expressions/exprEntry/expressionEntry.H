/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 Bernhard Gschaider
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::exprTools::expressionEntry

Description
    Convert dictionary entry to a stringified expression.

    The general OpenFOAM dictionary expansions will result in
    space-separated values. For example,

    \verbatim
        origin (0.21 0 0.01);

        condition  "mag(pos() - $centre) < 0.5";
    \endverbatim

    this will expand to the following:

    \verbatim
        condition  "mag(pos() - (0.21 0 0.01)) < 0.5";
    \endverbatim

    For these type of expressions, we'd would like better control.
    Using instead the special expansions, we can add an effective
    type cast.

    \verbatim
        condition  "mag(pos() - $[(vector)centre]) < 0.5";
    \endverbatim

    which will expand to the following:

    \verbatim
        condition  "mag(pos() - vector(0.21,0,0.01)) < 0.5";
    \endverbatim

SourceFiles
    expressionEntry.C
    expressionEntryI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_exprTools_expressionEntry_H
#define Foam_exprTools_expressionEntry_H

#include "exprString.H"
#include "primitiveEntry.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

/*---------------------------------------------------------------------------*\
                       Class expressionEntry Declaration
\*---------------------------------------------------------------------------*/

class expressionEntry
{
protected:

    //- Stringified version of data with comma-separated components.
    //- Uses prefix corresponding to the pTraits of Type.
    template<class Type>
    static string toExprStr(const Type& data);

    //- Comma-separated stringified version of primitiveEntry of Type.
    //- Prefix corresponding to the pTraits of Type
    template<class Type>
    static string toExprStr(ITstream& is);


public:

    //- Runtime type information
    TypeNameNoDebug("expressionEntry");

    declareRunTimeSelectionTable
    (
        autoPtr,
        expressionEntry,
        empty,
        (),
        ()
    );


    // Constructors

        //- Default construct
        expressionEntry() noexcept = default;


    // Selectors

        //- Return an entry to expression converter
        static autoPtr<expressionEntry> New(const word& name);


    //- Destructor
    virtual ~expressionEntry() = default;


    // Static Member Functions

        //- Generic concatenate tokens to space-separated string.
        inline static string evaluate(const entry& e);

        //- Inplace expand expression with dictionary variables/entries
        //
        // \par Expansion behaviour
        // - alternatives = True
        // - environment = True
        // - allow empty = True
        // - subDict = False
        // .
        static void inplaceExpand
        (
            std::string& s,
            const dictionary& dict
        );

        //- Expand expression with dictionary entries
        static expressions::exprString expand
        (
            const std::string& str,
            const dictionary& dict
        );


    // Member Functions

        //- To string. Normally with comma separators.
        virtual string toExpr(const entry& e) const
        {
            return evaluate(e);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace exprTools
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "expressionEntryI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
