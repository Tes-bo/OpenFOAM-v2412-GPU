/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

#include "fvExprDriverWriter.H"
#include "fvExprDriver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

defineTypeName(fvExprDriverWriter);

} // namespace expressions
} // namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::fvExprDriverWriter::fvExprDriverWriter
(
    const word& name,
    fvExprDriver& driver
)
:
    regIOobject
    (
        IOobject
        (
            name,
            driver.mesh().time().timeName(),
            "expressions",
            driver.mesh().time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    driver_(driver)
{
    if (headerOk())
    {
        readData(readStream("fvExprDriverWriter", true));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::fvExprDriverWriter::readData(Istream& is)
{
    dictionary dict(is);

    // driver_.readDict(is);

    driver_.getData(dict);

    return !is.bad();
}


bool Foam::expressions::fvExprDriverWriter::writeData(Ostream& os) const
{
    // driver_.writeDict(os);

    dictionary dict;
    driver_.prepareData(dict);
    dict.write(os, false);

    return os.good();
}


// ************************************************************************* //
