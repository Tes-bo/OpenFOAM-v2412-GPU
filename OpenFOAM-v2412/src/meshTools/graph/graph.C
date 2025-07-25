/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "graph.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "OSspecific.H"
#include "SubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef graph::writer graphWriter;
    defineTypeNameAndDebug(graphWriter, 0);
    defineRunTimeSelectionTable(graphWriter, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::graph::wordify(const Foam::string& sname)
{
    string wname = sname;
    wname.replace(" ", "_");
    wname.replace("(", "_");
    wname.replace(")", "");

    return word(wname);
}


void Foam::graph::readCurves(Istream& is)
{
    List<xy> xyData(is);

    x_.setSize(xyData.size());
    scalarField y(xyData.size());

    forAll(xyData, i)
    {
        x_[i] = xyData[i].x_;
        y[i] = xyData[i].y_;
    }

    set
    (
        wordify(yName_),
        new curve(wordify(yName_), curve::curveStyle::CONTINUOUS, y)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::graph::graph
(
    const string& title,
    const string& xName,
    const string& yName,
    const scalarField& x
)
:
    title_(title),
    xName_(xName),
    yName_(yName),
    x_(x)
{}


Foam::graph::graph
(
    const string& title,
    const string& xName,
    const string& yName,
    const scalarField& x,
    const scalarField& y
)
:
    title_(title),
    xName_(xName),
    yName_(yName),
    x_(x)
{
    set
    (
        wordify(yName),
        new curve(yName, curve::curveStyle::CONTINUOUS, y)
    );
}


Foam::graph::graph
(
    const string& title,
    const string& xName,
    const string& yName,
    Istream& is
)
:
    title_(title),
    xName_(xName),
    yName_(yName)
{
    readCurves(is);
}


Foam::graph::graph(Istream& is)
:
    title_(is),
    xName_(is),
    yName_(is)
{
    readCurves(is);
}


const Foam::scalarField& Foam::graph::y() const
{
    if (size() != 1)
    {
        FatalErrorInFunction
            << "y field requested for graph containing " << size()
            << "ys" << exit(FatalError);
    }

    return *begin()();
}


Foam::scalarField& Foam::graph::y()
{
    if (size() != 1)
    {
        FatalErrorInFunction
            << "y field requested for graph containing " << size()
            << "ys" << exit(FatalError);
    }

    return *begin()();
}


void Foam::graph::setXRange(const scalar x0, const scalar x1)
{
    if (x1 < x0)
    {
        FatalErrorInFunction
            << "When setting limits, x1 must be greater than x0" << nl
            << "    x0: " << x0 << nl
            << "    x1: " << x1 << nl
            << abort(FatalError);
    }

    label i0 = 0;
    label i1 = 0;

    forAll(x_, i)
    {
        if (x_[i] < x0)
        {
            i0 = i + 1;
        }
        if (x_[i] < x1)
        {
            i1 = i;
        }
    }

    label nX = i1 - i0 + 1;
    scalarField xNew(SubField<scalar>(x_, nX, i0));
    x_.transfer(xNew);

    forAllIters(*this, iter)
    {
        curve* c = iter();
        scalarField cNew(SubField<scalar>(*c, nX, i0));
        c->transfer(cNew);
    }
}


Foam::autoPtr<Foam::graph::writer> Foam::graph::writer::New
(
    const word& graphFormat
)
{
    if (!wordConstructorTablePtr_)
    {
        FatalErrorInFunction
            << "Graph writer table is empty"
            << exit(FatalError);
    }

    auto* ctorPtr = wordConstructorTable(graphFormat);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "graph",
            graphFormat,
            *wordConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<graph::writer>(ctorPtr());
}


void Foam::graph::writer::writeXY
(
    const scalarField& x,
    const scalarField& y,
    Ostream& os
) const
{
    forAll(x, xi)
    {
        os << setw(10) << x[xi] << token::SPACE << setw(10) << y[xi]<< endl;
    }
}


void Foam::graph::writeTable(Ostream& os) const
{
    forAll(x_, xi)
    {
        os  << setw(10) << x_[xi];

        forAllConstIters(*this, iter)
        {
            os  << token::SPACE << setw(10) << (*iter())[xi];
        }
        os  << endl;
    }
}


void Foam::graph::write(Ostream& os, const word& format) const
{
    writer::New(format)().write(*this, os);
}


void Foam::graph::write(const fileName& pName, const word& format) const
{
    autoPtr<writer> writer(writer::New(format));

    OFstream graphFile(pName + '.' + writer().ext());

    if (graphFile.good())
    {
        write(graphFile, format);
    }
    else
    {
        WarningInFunction
            << "Could not open graph file " << graphFile.name()
            << endl;
    }
}


void Foam::graph::write
(
    const fileName& path,
    const word& name,
    const word& format
) const
{
    mkDir(path);
    write(path/name, format);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const graph& g)
{
    g.writeTable(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
