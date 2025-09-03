/*---------------------------------------------------------------------------*\
  =========                 |
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

\*---------------------------------------------------------------------------*/

#include "graphFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "labelVector.H"
#include "FlatOutput.H"
#include "SVGTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(graphFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        graphFunctionObject,
        dictionary
    );
}
}

// 'Muted' colour scheme from https://personal.sron.nl/~pault/ (12.07.24)
Foam::wordList Foam::functionObjects::graphFunctionObject::defaultColours
({
    "#CC6677",
    "#332288",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499"
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::graphFunctionObject::getValue
(
    const label objecti,
    label& valuei
)
{
    const word& object = objects_[objecti];
    const word& entry = entries_[objecti];

    Type result;
    if (!this->getObjectResult(object, entry, result))
    {
        return false;
    }

    auto& cols = objectToCol_[objecti];
    if (cols.empty())
    {
        for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
        {
            cols.push_back(valuei++);
            values_.push_back(DynamicList<scalar>());
        }
    }

    for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
    {
        scalar v = component(result, d);

        if (logScaleY_)
        {
            v = (v < SMALL) ? 1 : log10(v);
        }

        values_[cols[d]].push_back(v);
    }

    return true;
}


Foam::label Foam::functionObjects::graphFunctionObject::setAxisProps
(
    const bool logScale,
    scalar& xmin,
    scalar& xmax,
    scalar& xtick
) const
{
    DebugInfo
        << "1 -- xmin:" << xmin << " xmax:" << xmax
        << " xtick:" << xtick << endl;

    /*
    Divisions Based on (12.07.24):
    https://peltiertech.com/calculate-nice-axis-scales-in-your-excel-worksheet
    */

    const scalar range = xmax - xmin;
    const scalar eps = 0.01*range;

    // Extend xmin and xmax by eps
    if (mag(xmin) < SMALL)
    {
        xmin = 0;
    }
    else
    {
        xmin = (xmin > 0) ? max(0, xmin - eps) : xmin - eps;
    }

    if (mag(xmax) < SMALL)
    {
        xmax = mag(xmin) < SMALL ? 1 : 0;
    }
    else
    {
        xmax = (xmax < 0) ? min(0, xmax + eps) : xmax + eps;
    }

    DebugInfo
        << "2 -- xmin:" << xmin << " xmax:" << xmax
        << " xtick:" << xtick << endl;

    auto lookup = [](const scalar x) -> scalar
    {
        if (x < 2.5) { return 0.2; }
        if (x < 5.0) { return 0.5; }
        if (x < 10.0) { return 2.0; }
        return 10.0;
    };

    const scalar power = log10(range);
    const scalar factor = pow(10, power - floor(power));

    xtick = lookup(factor)*pow(10, floor(power));
    xmin = xtick*floor(xmin/xtick);
    xmax = xtick*(floor(xmax/xtick) + 1);

    // Convert ticks to integer powers of 10 for log scales
    if (logScale)
    {
        xmin = floor(xmin);
        xmax = ceil(xmax);
        xtick = 1;
    }

    DebugInfo
        << "power:" << power << " factor:" << factor
        << " xmin:" << xmin << " xmax:" << xmax
        << " xtick:" << xtick << endl;

    return round((xmax - xmin)/xtick);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::graphFunctionObject::graphFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    writeFile(runTime, name, typeName, dict, true, ".svg"),
    objects_(),
    entries_(),
    titles_(),
    colours_(),
    dashes_(),
    times_(),
    values_(),
    objectToCol_(),
    xMin_(dict.getOrDefault<scalar>("xMin", GREAT)),
    xMax_(dict.getOrDefault<scalar>("xMax", GREAT)),
    yMin_(dict.getOrDefault<scalar>("yMin", GREAT)),
    yMax_(dict.getOrDefault<scalar>("yMax", GREAT)),
    xlabel_(dict.getOrDefault<string>("xlabel", "Iteration/Time")),
    ylabel_(dict.getOrDefault<string>("ylabel", "Property")),
    width_(dict.getOrDefault<label>("width", 800)),
    height_(dict.getOrDefault<label>("height", 600)),
    strokeWidth_(dict.getOrDefault<label>("strokeWidth", 2)),
    logScaleX_(dict.getOrDefault<bool>("logScaleX", false)),
    logScaleY_(dict.getOrDefault<bool>("logScaleY", false)),
    drawGrid_(dict.getOrDefault<bool>("drawGrid", true))
{
    const dictionary& functions = dict.subDict("functions");
    objects_.setSize(functions.size());
    entries_.setSize(functions.size());
    titles_.setSize(functions.size());
    colours_.setSize(functions.size());
    dashes_.setSize(functions.size());
    objectToCol_.setSize(functions.size());

    label defaultColouri = 0;
    label entryi = 0;

    for (const auto& e : functions)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(functions)
                << "Functions must be provided in dictionary format"
                << exit(FatalIOError);
        }

        const dictionary& d = e.dict();
        objects_[entryi] = d.get<word>("object");
        entries_[entryi] = d.get<word>("entry");
        titles_[entryi] = d.getOrDefault<string>("title", e.keyword());

        labelVector colour;
        if (d.readIfPresent("colour", colour))
        {
            // Warn/error if outside 0-255 range?
            colour[0] = min(255, max(0, colour[0]));
            colour[1] = min(255, max(0, colour[1]));
            colour[2] = min(255, max(0, colour[2]));

            OStringStream oss;
            oss << "rgb" << flatOutput(colour, FlatOutput::ParenComma{});
            colours_[entryi] = oss.str();
        }
        else
        {
            colours_[entryi] = defaultColours[defaultColouri++];
            if (defaultColouri == defaultColours.size())
            {
                // Lots of lines to plot - exhausted list of default colours.
                // Restarting ...
                defaultColouri = 0;
            }
        }

        {
            labelList dashes;
            if (d.readIfPresent("dashes", dashes))
            {
                OStringStream oss;
                oss << flatOutput(dashes, FlatOutput::BareSpace{});
                dashes_[entryi] = oss.str();
            }
            else
            {
                // Solid line
                dashes_[entryi] = "0";
            }
        }

        ++entryi;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::graphFunctionObject::execute()
{
    if (!Pstream::master()) return true;

    scalar& graphTime = times_.emplace_back(time_.timeOutputValue());

    if (logScaleX_)
    {
        graphTime = log10(max(graphTime, SMALL));
    }

    label valuei = 0;
    forAll(objects_, objecti)
    {
        bool ok =
            getValue<label>(objecti, valuei)
         || getValue<scalar>(objecti, valuei)
         || getValue<vector>(objecti, valuei)
         || getValue<sphericalTensor>(objecti, valuei)
         || getValue<symmTensor>(objecti, valuei)
         || getValue<tensor>(objecti, valuei);

        if (!ok)
        {
            // Entry not found
            Log << type() << " " << name() << " execute: "
                << "Unable to get value for object:" << objects_[objecti]
                << " entry:" << entries_[objecti] << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::graphFunctionObject::write()
{
    if (!Pstream::master()) return true;
    // DebugVar(values_);

    auto filePtr = newFileAtTime(name(), time().value());
    auto& os = filePtr();

    scalar ymin = GREAT;
    scalar ymax = -GREAT;

    bool valid = false;
    for (const auto& data : values_)
    {
        for (const auto& value : data)
        {
            ymin = min(ymin, value);
            ymax = max(ymax, value);
            valid = true;
        }
    }

    // Early exit if there is no data
    if (!valid)
    {
        Log << type() << " " << name() << " write:" << nl
            << "    No data to plot - skipping" << nl;

        // Empty graph
        os  << SVG::header(width_, height_) << nl << SVG::end << endl;

        return false;
    }

    auto applyLimits = [](const scalar val, const scalar lim, const bool lg)
    {
        if (lim < 0.99*GREAT)
        {
            return lg ? log10(lim) : lim;
        }

        return val;
    };


    // Set y axis limits if user-supplied
    ymin = applyLimits(ymin, yMin_, logScaleY_);
    ymax = applyLimits(ymax, yMax_, logScaleY_);


    scalar ytick = 0;
    const label ny = setAxisProps(logScaleY_, ymin, ymax, ytick);

    const scalar border = 0.1;
    const scalar w = width_*(1.0 - 2*border);
    const scalar h = height_*(1.0 - 2*border);

    // Set x axis limits if user-supplied
    scalar xmin = applyLimits(0, xMin_, logScaleX_);
    scalar xmax = applyLimits(max(times_), xMax_, logScaleX_);

    // Set x axis properties; return the number of tic divisions
    scalar xtick = 0;
    const label nx = setAxisProps(logScaleX_, xmin, xmax, xtick);

    // Top pixel co-ordinate
    auto top = [=](const scalar y)
    {
        const scalar ratio = (y - ymin)/(ymax - ymin + ROOTVSMALL);
        return round(height_ - ratio*h - border*height_);
    };

    // Left pixel co-ordinate
    auto left = [=](const scalar x)
    {
        return round(x/(xmax - xmin + ROOTVSMALL)*w + border*width_);
    };

    const scalar fontpx = min(20, h/(2*values_.size()));
    const scalar fontdy = 1.5*fontpx;

    // Legend - top right: text (right aligned), coloured line (fixed positions)
    const label legendLineRight = border*width_ + w - fontpx;
    const label legendLineLeft = legendLineRight - 0.5*border*width_;
    const label legendLabelRight = legendLineLeft - 0.5*fontpx;

    // Graph box and tick colour
    const word colour = "rgb(105,105,105)";

    os  << SVG::header(width_, height_) << nl;

    // Graph bounding box
    SVG::element bounds("rect", {{"fill", "none"}, {"stroke", colour}});
    bounds.addAttr("x", round(border*width_));
    bounds.addAttr("y", round(border*height_));
    bounds.addAttr("width", round(w));
    bounds.addAttr("height", round(h));
    os  << bounds << bounds.end << nl;


    // X axis label
    os  << SVG::text
        (
            xlabel_,
            0.5*width_,
            height_ - 0.5*(border*height_) + fontpx,
            {{"font-size", Foam::name(1.2*fontpx)}},
            "middle"
        )
        << nl;

    // Y axis label - text rotated
    SVG::text ytext
    (
        ylabel_,
        0,
        0,
        {{"font-size", Foam::name(1.2*fontpx)}},
        "middle"
    );
    ytext.addAttr("alignment-baseline", "middle");
    ytext.addAttrStr
    (
        "transform",
        "translate(" + Foam::name(left(xmin) - 3*fontpx) + ","
      + Foam::name(0.5*height_) + ") rotate(270)"
    );
    os  << ytext << nl;

    const label dTick = 0.2*fontpx;

    // Background grid
    if (drawGrid_)
    {
        const word colourGrid = "rgb(200,200,200)";

        for (label i = 1; i < nx; ++i)
        {
            const label x = left(xmin + i*xtick);
            const label y1 = top(ymin);
            const label y2 = top(ymax);

            // Dashed grid lines
            os  << SVG::line
                (
                    x,
                    y1,
                    x,
                    y2,
                    {
                        {"stroke", colourGrid},
                        {"stroke-width", "1"},
                        {"stroke-dasharray", "4"}
                    }
                ) << nl;
        }

        for (label i = 1; i < ny; ++i)
        {
            const label y = top(ymin + i*ytick);
            const label x1 = left(xmin);
            const label x2 = left(xmax);

            // Dashed grid lines
            os  << SVG::line
                (
                    x1,
                    y,
                    x2,
                    y,
                    {
                        {"stroke", colourGrid},
                        {"stroke-width", "1"},
                        {"stroke-dasharray", "4"}
                    }
                ) << nl;
        }
    }

    // Axis labels
    for (label i = 0; i <= nx; ++i)
    {
        const scalar v = xmin + i*xtick;
        const label x = left(v);
        const scalar y0 = ymin;
        const label y1 = top(y0);
        const label y2 = y1 + dTick;
        const string tickLabel = logScaleX_
           ? "<tspan>10<tspan style=\"font-size:"
           + Foam::name(label(0.75*fontpx))
           + "px\" dy=\"" + Foam::name(-0.4*fontpx) + "\">"
           + Foam::name(v)
           + "</tspan></tspan>"
           : Foam::name(v);

        // Ticks
        os  << SVG::line
            (
                x,
                y1,
                x,
                y2,
                {
                    {"stroke", colour},
                    {"stroke-width", Foam::name(strokeWidth_)}
                }
            ) << nl;

        // Labels
        os  << SVG::text
            (
                tickLabel,
                x,
                y2 + 1.25*fontpx,
                {{"font-size", Foam::name(fontpx)}},
                "middle"
            ) << nl;
    }
    for (label i = 0; i <= ny; ++i)
    {
        const scalar v = ymin + i*ytick;
        const label y = top(v);
        const label y2 = y + 0.4*fontpx;
        const scalar x0 = xmin;
        const label x1 = left(x0);
        const label x2 = x1 - dTick;
        const string tickLabel = logScaleY_
          ? "<tspan>10<tspan style=\"font-size:"
          + Foam::name(label(0.6*fontpx))
          + "px\" dy=\"" + Foam::name(-0.4*fontpx) + "\">"
          + Foam::name(v)
          + "</tspan></tspan>"
          : Foam::name(v);

        // Ticks
        os  << SVG::line
            (
                x1,
                y,
                x2,
                y,
                {{"stroke", colour},{"stroke-width",  "1"}}
            ) << nl;

        // Labels
        os  << SVG::text
            (
                tickLabel,
                x2 - 0.5*fontpx,
                y2,
                {{"font-size", Foam::name(fontpx)}},
                "end"
            ) << nl;
    }


    forAll(objects_, objecti)
    {
        const word& colour = colours_[objecti];

        const auto& cols = objectToCol_[objecti];
        for (const label c : cols)
        {
            const word cmpt = cols.size() > 1 ? Foam::name(c) : "";

            label legendTop = border*height_ + fontdy*(c+1);
            os  << SVG::text
                (
                    titles_[objecti] + cmpt,
                    legendLabelRight,
                    legendTop,
                    {{"font-size", Foam::name(fontpx)}},
                    "end"
                ) << nl;

            os  << SVG::line
                (
                    legendLineLeft,
                    legendTop - 0.5*fontpx,
                    legendLineRight,
                    legendTop - 0.5*fontpx,
                    {{"stroke", colour},{"stroke-width", "2"}},
                    {{"stroke-dasharray", dashes_[objecti]}}
                ) << nl;


            os  << "<path d=\"";
            const auto& data = values_[c];
            bool firstPoint = true;
            forAll(data, i)
            {
                const scalar t = times_[i];
                const scalar v = data[i];

                if ((v > ymin) && (v < ymax))
                {
                    if (firstPoint)
                    {
                        os  << " M ";
                    }
                    else
                    {
                        os  << " L ";
                    }

                    os  << left(t) << ' ' << top(v);

                    firstPoint = false;
                }
                else
                {
                    firstPoint = true;
                }
            }

            os  << "\""
                << " style=\"stroke:" << colour << ";"
                << " fill:none; stroke-width:2;\""
                << " stroke-dasharray=\"" << dashes_[objecti].c_str() << "\" />"
                << nl;
        }
    }

    os  << SVG::end << endl;

    Log << type() << " " << name() <<  " write:" << nl
        << "    Written file " << os.name() << nl << endl;

    return true;
}


// ************************************************************************* //