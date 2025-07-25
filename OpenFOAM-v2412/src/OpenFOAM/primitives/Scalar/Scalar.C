/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const pTraits<Scalar>::typeName = "scalar";
const char* const pTraits<Scalar>::componentNames[] = { "" };

const Scalar pTraits<Scalar>::zero = 0.0;
const Scalar pTraits<Scalar>::one = 1.0;
const Scalar pTraits<Scalar>::min = -ScalarVGREAT;
const Scalar pTraits<Scalar>::max = ScalarVGREAT;
const Scalar pTraits<Scalar>::rootMin = -ScalarROOTVGREAT;
const Scalar pTraits<Scalar>::rootMax = ScalarROOTVGREAT;
const Scalar pTraits<Scalar>::vsmall = ScalarVSMALL;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pTraits<Scalar>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * IO/Conversion * * * * * * * * * * * * * * * //

word name(const Scalar val)
{
    // Caution std::to_string(double) is locale sensitive!
    std::ostringstream buf;
    buf << val;

    return word(buf.str(), false);  // Needs no stripping
}


Scalar ScalarRead(const char* buf)
{
    char* endptr = nullptr;
    errno = 0;
    const auto parsed = ScalarConvert(buf, &endptr);

    const parsing::errorType err =
    (
        (parsed < -ScalarVGREAT || parsed > ScalarVGREAT)
      ? parsing::errorType::RANGE
      : parsing::checkConversion(buf, endptr)
    );

    if (err != parsing::errorType::NONE)
    {
        FatalIOErrorInFunction("unknown")
            << parsing::errorNames[err] << " '" << buf << "'"
            << exit(FatalIOError);
    }

    // Round underflow to zero
    return
    (
        (parsed > -ScalarVSMALL && parsed < ScalarVSMALL)
      ? 0
      : Scalar(parsed)
    );
}


bool ScalarRead(const char* buf, Scalar& val)
{
    char* endptr = nullptr;
    errno = 0;
    const auto parsed = ScalarConvert(buf, &endptr);

    // Round underflow to zero
    val =
    (
        (parsed >= -ScalarVSMALL && parsed <= ScalarVSMALL)
      ? 0
      : Scalar(parsed)
    );

    return
    (
        (parsed < -ScalarVGREAT || parsed > ScalarVGREAT)
      ? false
      : (parsing::checkConversion(buf, endptr) == parsing::errorType::NONE)
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Scalar ScalarRead(Istream& is)
{
    Scalar val(0);
    is  >> val;

    return val;
}


Istream& operator>>(Istream& is, Scalar& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get scalar value"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    // Accept separated '-' (or '+') while expecting a number.
    // This can arise during dictionary expansions (Eg, -$value)

    char prefix = 0;
    if (t.isPunctuation())
    {
        prefix = t.pToken();
        if (prefix == token::PLUS || prefix == token::MINUS)
        {
            is >> t;
        }
    }

    if (t.isNumber())
    {
        val =
        (
            (prefix == token::MINUS)
          ? (0 - t.number())
          : t.number()
        );
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected scalar value, found ";
        if (prefix == token::PLUS || prefix == token::MINUS)
        {
            FatalIOError << '\'' << prefix << "' followed by ";
        }
        FatalIOError << t.info() << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Ostream& operator<<(Ostream& os, const Scalar val)
{
    os.write(val);
    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
