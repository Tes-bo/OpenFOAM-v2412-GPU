/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

InNamespace
    Foam::stringOps

Description
    Specialized string sorting.

SourceFiles
    stringOpsSort.C

\*---------------------------------------------------------------------------*/

#ifndef stringOpsSort_H
#define stringOpsSort_H

#include "stringOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stringOps
{

    //- 'Natural' compare for C-strings
    //  Uses algorithm and code from Jan-Marten Spit <jmspit@euronet.nl>
    //
    //  In the 'natural' comparison, strings are compared alphabetically
    //  and numerically. Thus 'file010.txt' sorts after 'file2.txt'
    //
    //  \param s1 left string
    //  \param s2 right string
    //  \return -1 when s1 < s2, 0 when s1 == s2, 1 when s1 > s2
    int natstrcmp(const char* s1, const char* s2);


    //- Encapsulation of natural order sorting for algorithms
    struct natural_sort
    {
        //- Natural compare for std::string
        //  \return -1 when s1 < s2, 0 when s1 == s2, 1 when s1 > s2
        static inline int compare
        (
            const std::string& s1,
            const std::string& s2
        )
        {
            return natstrcmp(s1.data(), s2.data());
        }

        //- Default (forward) natural sorting
        bool operator()(const std::string& s1, const std::string& s2) const
        {
            return natural_sort::compare(s1, s2) < 0;
        }

        //- Reverse natural sorting
        struct reverse
        {
            //- Reverse natural sorting
            bool operator()(const std::string& s1, const std::string& s2) const
            {
                return natural_sort::compare(s1, s2) > 0;
            }
        };


        //- A list compare binary predicate for natural sort
        template<class T>
        struct less
        {
            const UList<T>& values;

            less(const UList<T>& list)
            :
                values(list)
            {}

            bool operator()(const label a, const label b) const
            {
                return natural_sort::compare(values[a], values[b]) < 0;
            }
        };


        //- A list compare binary predicate for reverse natural sort
        template<class T>
        struct greater
        {
            const UList<T>& values;

            greater(const UList<T>& list)
            :
                values(list)
            {}

            bool operator()(const label a, const label b) const
            {
                return natural_sort::compare(values[a], values[b]) > 0;
            }
        };
    };

} // End namespace stringOps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
