/*---------------------------------------------------------------------------*\
  =========                |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::SVG

Description
    Collection of tools to generate SVG strings

SourceFiles
    SVGTools.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SVGTools_H
#define Foam_SVGTools_H

#include "Ostream.H"
#include "OStringStream.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Namespace SVG Declaration
\*---------------------------------------------------------------------------*/

namespace SVG
{
    typedef std::pair<const char*, string> entryType;

    struct element;
    Ostream& operator<<(Ostream& os, const element& e);

    // Base SVG element
    struct element
    {
        const word key_;
        DynamicList<entryType> styles_;
        DynamicList<entryType> elems_;

        element
        (
            const word& key,
            const std::initializer_list<entryType>& styles  = {},
            const std::initializer_list<entryType>& elems = {}
        )
        :
            key_(key),
            styles_(styles),
            elems_(elems)
        {}

        template<class Type>
        void addAttr(const char* key, const Type& value)
        {
            OStringStream oss;
            oss << value;
            elems_.push_back(entryType(key, oss.str().c_str()));
        }

        void addAttrStr(const char* key, const string& str)
        {
            elems_.push_back(entryType(key, str.c_str()));
        }

        friend Ostream& operator<<(Ostream& os, const element& ele)
        {
            os  << "<" << ele.key_;

            for (const auto& e : ele.elems_)
            {
                os << " " << e.first << "=" << e.second;
            }

            os  << " style=\"";
            for (const auto& s : ele.styles_)
            {
                os << s.first << ":" << s.second.c_str() << ";";
            }

            os << "\">";

            return os;
        }

        const word end = "</" + key_ + ">";
    };


    struct text;
    Ostream& operator<<(Ostream& os, const text& t);

    // Text
    struct text
    :
        element
    {
        const string text_;

        text
        (
            const string text,
            const label left,
            const label top,
            const std::initializer_list<entryType>& styles = {},
            const word anchor = "middle",
            const std::initializer_list<entryType>& elems = {}
        )
        :
            element("text", styles, elems),
            text_(text)
        {
            elems_.push_back(entryType("x", Foam::name(left)));
            elems_.push_back(entryType("y", Foam::name(top)));
            elems_.push_back(entryType("text-anchor", anchor));
            elems_.push_back
            (
                entryType("font-family", "Arial, Helvetica, sans-serif")
            );
        }


        friend Ostream& operator<<(Ostream& os, const text& t)
        {
            // element::operator<<(os, t);
            os  << static_cast<const element&>(t);

            os  << t.text_.c_str();

            os  << t.end;

            return os;
        }
    };


    struct line;
    Ostream& operator<<(Ostream& os, const line& l);

    // Line
    struct line
    :
        element
    {
        line
        (
            const label x1,
            const label y1,
            const label x2,
            const label y2,
            const std::initializer_list<entryType>& styles = {},
            const std::initializer_list<entryType>& elems = {}
        )
        :
            element("line", styles, elems)
        {
            elems_.push_back(entryType("x1", Foam::name(x1)));
            elems_.push_back(entryType("y1", Foam::name(y1)));
            elems_.push_back(entryType("x2", Foam::name(x2)));
            elems_.push_back(entryType("y2", Foam::name(y2)));
        }


        friend Ostream& operator<<(Ostream& os, const line& l)
        {
            // element::operator<<(os, l);
            os  << static_cast<const element&>(l);
            os  << l.end;

            return os;
        }
    };

    struct header;
    Ostream& operator<<(Ostream& os, const header& h);

    // Header
    struct header
    {
        label width_;
        label height_;

        header(const label width, const label height)
        :
            width_(width),
            height_(height)
        {}

        friend Ostream& operator<<(Ostream& os, const header& h)
        {
            os  << "<svg viewBox=\"0 0 " << h.width_ << ' ' << h.height_ << "\""
                << " xmlns=\"http://www.w3.org/2000/svg\""
                << " xmlns:xlink=\"http://www.w3.org/1999/xlink\""
                << " xmlns:bx=\"https://www.boxy-svg.com/bx\">";

            return os;
        }
    };

    // Close SVG element
    const char* end  = "</svg>";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SVG
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //