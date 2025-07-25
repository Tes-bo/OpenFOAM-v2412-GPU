/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "dimensionSet.H"
#include "dimensionedScalar.H"
#include "IOstreams.H"
#include <limits>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dimensionSet::dimensionSet
(
    const word& entryName,
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
:
    exponents_(Zero)
{
    readEntry(entryName, dict, readOpt);
}


Foam::dimensionSet::dimensionSet(Istream& is)
{
    is >> *this;
}


Foam::dimensionSet::tokeniser::tokeniser(Istream& is)
:
    is_(is),
    tokens_(100),
    start_(0),
    size_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dimensionSet::tokeniser::push(const token& t)
{
    const label end = (start_+size_)%tokens_.size();
    tokens_[end] = t;
    if (size_ == tokens_.size())
    {
        start_ = tokens_.fcIndex(start_);
    }
    else
    {
        ++size_;
    }
}


Foam::token Foam::dimensionSet::tokeniser::pop()
{
    token t = tokens_[start_];
    start_ = tokens_.fcIndex(start_);
    --size_;
    return t;
}


void Foam::dimensionSet::tokeniser::unpop(const token& t)
{
    ++size_;
    start_ = tokens_.rcIndex(start_);
    tokens_[start_] = t;
}


bool Foam::dimensionSet::tokeniser::hasToken() const
{
    return size_ || is_.good();
}


bool Foam::dimensionSet::tokeniser::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // div
     && c != ';'   // end statement
     && c != '{'   // beg subdict
     && c != '}'   // end subdict
     && c != '('   // beg expr
     && c != ')'   // end expr
     && c != '['   // beg dim
     && c != ']'   // end dim
     && c != '^'   // power
     && c != '*'   // mult
    );
}


Foam::label Foam::dimensionSet::tokeniser::priority(const token& t)
{
    if (t.isPunctuation())
    {
        if
        (
            t.pToken() == token::MULTIPLY
         || t.pToken() == token::DIVIDE
        )
        {
            return 2;
        }
        else if (t.pToken() == '^')
        {
            return 3;
        }
    }

    // Default priority
    return 0;
}


void Foam::dimensionSet::tokeniser::splitWord(const word& w)
{
    size_t start = 0;
    for (size_t i=0; i<w.size(); ++i)
    {
        if (!valid(w[i]))
        {
            if (i > start)
            {
                const word subWord = w.substr(start, i-start);
                if (isdigit(subWord[0]) || subWord[0] == token::SUBTRACT)
                {
                    push(token(readScalar(subWord)));
                }
                else
                {
                    push(token(subWord));
                }
            }
            if (w[i] != token::SPACE)
            {
                if (isdigit(w[i]))
                {
                    // Single digit: as scalar value
                    const scalar val = (w[i] - '0');
                    push(token(val));
                }
                else
                {
                    push(token(token::punctuationToken(w[i])));
                }
            }
            start = i+1;
        }
    }
    if (start < w.size())
    {
        const word subWord = w.substr(start);
        if (isdigit(subWord[0]) || subWord[0] == token::SUBTRACT)
        {
            push(token(readScalar(subWord)));
        }
        else
        {
            push(token(subWord));
        }
    }
}


Foam::token Foam::dimensionSet::tokeniser::nextToken()
{
    if (size_ == 0)
    {
        token t(is_);
        if (t.isWord())
        {
            splitWord(t.wordToken());
            return pop();
        }
        else
        {
            return t;
        }
    }
    else
    {
        return pop();
    }
}


void Foam::dimensionSet::tokeniser::putBack(const token& t)
{
    if (size_ == 0)
    {
        push(t);
    }
    else
    {
        unpop(t);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dimensionSet::round(const scalar tol)
{
    scalar integralPart;
    for (scalar& val : exponents_)
    {
        const scalar fractionalPart = std::modf(val, &integralPart);

        if (mag(fractionalPart-1.0) <= tol)
        {
            val = 1.0+integralPart;
        }
        else if (mag(fractionalPart+1.0) <= tol)
        {
            val = -1.0+integralPart;
        }
        else if (mag(fractionalPart) <= tol)
        {
            val = integralPart;
        }
    }
}


Foam::dimensionedScalar Foam::dimensionSet::parse
(
    const label lastPrior,
    tokeniser& tis,
    const HashTable<dimensionedScalar>& readSet
) const
{
    dimensionedScalar ds("", dimless, 1);

    // Get initial token
    token nextToken(tis.nextToken());

    // Store type of last token read. Used to detect two consecutive
    // symbols and assume multiplication
    bool haveReadSymbol = false;


    while (true)
    {
        if (nextToken.isWord())
        {
            const word& unitName = nextToken.wordToken();
            const dimensionedScalar& unitDim = readSet[unitName];
            ds.dimensions() *= unitDim.dimensions();
            ds.value() *= unitDim.value();
            haveReadSymbol = true;
        }
        else if (nextToken.isNumber())
        {
            // no dimensions, just value
            ds.value() *= nextToken.number();
            haveReadSymbol = true;
        }
        else if (nextToken.isPunctuation())
        {
            label nextPrior = tokeniser::priority(nextToken);

            if (nextToken.pToken() == token::BEGIN_SQR)
            {
                // No idea when this will happen
                tis.putBack(nextToken);
                return ds;
            }
            else if (nextToken.pToken() == token::END_SQR)
            {
                tis.putBack(nextToken);
                return ds;
            }
            else if (nextToken.pToken() == token::BEGIN_LIST)
            {
                dimensionedScalar sub(parse(nextPrior, tis, readSet));

                token t = tis.nextToken();
                if (!t.isPunctuation()  || t.pToken() !=  token::END_LIST)
                {
                    FatalIOErrorInFunction(tis.stream())
                        << "Illegal token " << t << exit(FatalIOError);
                }

                ds.dimensions() *= sub.dimensions();
                ds.value() *= sub.value();

                haveReadSymbol = true;
            }
            else if (nextToken.pToken() == token::END_LIST)
            {
                tis.putBack(nextToken);
                return ds;
            }
            else if (nextToken.pToken() == token::MULTIPLY)
            {
                if (nextPrior > lastPrior)
                {
                    dimensionedScalar sub(parse(nextPrior, tis, readSet));

                    ds.dimensions() *= sub.dimensions();
                    ds.value() *= sub.value();
                }
                else
                {
                    // Restore token
                    tis.putBack(nextToken);
                    return ds;
                }
                haveReadSymbol = false;
            }
            else if (nextToken.pToken() == token::DIVIDE)
            {
                if (nextPrior > lastPrior)
                {
                    dimensionedScalar sub(parse(nextPrior, tis, readSet));

                    ds.dimensions() /= sub.dimensions();
                    ds.value() /= sub.value();
                }
                else
                {
                    tis.putBack(nextToken);
                    return ds;
                }
                haveReadSymbol = false;
            }
            else if (nextToken.pToken() == '^')
            {
                if (nextPrior > lastPrior)
                {
                    dimensionedScalar expon(parse(nextPrior, tis, readSet));

                    ds.dimensions().reset(pow(ds.dimensions(), expon.value()));
                    // Round to nearest integer if close to it
                    ds.dimensions().round(10*smallExponent);
                    ds.value() = Foam::pow(ds.value(), expon.value());
                }
                else
                {
                    tis.putBack(nextToken);
                    return ds;
                }
                haveReadSymbol = false;
            }
            else
            {
                FatalIOErrorInFunction(tis.stream())
                    << "Illegal token " << nextToken << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction(tis.stream())
                << "Illegal token " << nextToken << exit(FatalIOError);
        }


        if (!tis.hasToken())
        {
            break;
        }

        nextToken = tis.nextToken();
        if (nextToken.error())
        {
            break;
        }

        if (haveReadSymbol && (nextToken.isWord() || nextToken.isNumber()))
        {
            // Two consecutive symbols. Assume multiplication
            tis.putBack(nextToken);
            nextToken = token(token::MULTIPLY);
        }
    }

    return ds;
}


bool Foam::dimensionSet::readEntry
(
    const word& entryName,
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (readOpt == IOobjectOption::NO_READ)
    {
        return false;
    }

    const entry* eptr = dict.findEntry(entryName, keyType::LITERAL);

    if (eptr)
    {
        const entry& e = *eptr;
        ITstream& is = e.stream();

        is >> *this;

        e.checkITstream(is);

        return true;
    }
    else if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(dict)
            << "Entry '" << entryName << "' not found in dictionary "
            << dict.relativeName() << nl
            << exit(FatalIOError);
    }

    return false;
}


Foam::Istream& Foam::dimensionSet::read
(
    Istream& is,
    scalar& multiplier,
    const HashTable<dimensionedScalar>& readSet
)
{
    multiplier = 1.0;

    // Read beginning of dimensionSet
    token startToken(is);

    if (startToken != token::BEGIN_SQR)
    {
        FatalIOErrorInFunction(is)
            << "Expected a '" << token::BEGIN_SQR << "' in dimensionSet\n"
            << "in stream " << is.info() << nl
            << exit(FatalIOError);
    }

    // Read next token
    token nextToken(is);

    if (!nextToken.isNumber())
    {
        is.putBack(nextToken);

        tokeniser tis(is);

        dimensionedScalar ds(parse(0, tis, readSet));

        multiplier = ds.value();
        exponents_ = ds.dimensions().values();
    }
    else
    {
        // Read first five dimensions
        exponents_[dimensionSet::MASS] = nextToken.number();
        for (int d=1; d < dimensionSet::CURRENT; ++d)
        {
            is >> exponents_[d];
        }

        // Read next token
        token nextToken(is);

        // If next token is another number
        // read last two dimensions
        // and then read another token for the end of the dimensionSet
        if (nextToken.isNumber())
        {
            exponents_[dimensionSet::CURRENT] = nextToken.number();
            is >> nextToken;
            exponents_[dimensionSet::LUMINOUS_INTENSITY] = nextToken.number();
            is >> nextToken;
        }
        else
        {
            exponents_[dimensionSet::CURRENT] = 0;
            exponents_[dimensionSet::LUMINOUS_INTENSITY] = 0;
        }

        // Check end of dimensionSet
        if (nextToken != token::END_SQR)
        {
            FatalIOErrorInFunction(is)
                << "Expected a '" << token::END_SQR << "' in dimensionSet\n"
                << "in stream " << is.info() << nl
                << exit(FatalIOError);
        }
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Istream& Foam::dimensionSet::read
(
    Istream& is,
    scalar& multiplier
)
{
    return read(is, multiplier, unitSet());
}


Foam::Istream& Foam::dimensionSet::read
(
    Istream& is,
    scalar& multiplier,
    const dictionary& readSet
)
{
    multiplier = 1.0;

    // Read beginning of dimensionSet
    token startToken(is);

    if (startToken != token::BEGIN_SQR)
    {
        FatalIOErrorInFunction(is)
            << "Expected a '" << token::BEGIN_SQR << "' in dimensionSet\n"
            << "in stream " << is.info() << nl
            << exit(FatalIOError);
    }

    // Read next token
    token nextToken(is);

    if (nextToken.isWord())
    {
        bool continueParsing = true;
        do
        {
            word symbolPow = nextToken.wordToken();
            if (symbolPow.back() == token::END_SQR)
            {
                symbolPow.resize(symbolPow.size()-1);
                continueParsing = false;
            }


            // Parse unit
            dimensionSet symbolSet;  // dimless

            const auto index = symbolPow.find('^');
            if (index != std::string::npos)
            {
                const word symbol = symbolPow.substr(0, index);
                const scalar exponent = readScalar(symbolPow.substr(index+1));

                dimensionedScalar s;
                s.read(readSet.lookup(symbol, keyType::LITERAL), readSet);

                symbolSet.reset(pow(s.dimensions(), exponent));

                // Round to nearest integer if close to it
                symbolSet.round(10*smallExponent);
                multiplier *= Foam::pow(s.value(), exponent);
            }
            else
            {
                dimensionedScalar s;
                s.read(readSet.lookup(symbolPow, keyType::LITERAL), readSet);

                symbolSet.reset(s.dimensions());
                multiplier *= s.value();
            }

            // Add dimensions without checking
            for (int i=0; i < dimensionSet::nDimensions; ++i)
            {
                exponents_[i] += symbolSet[i];
            }

            if (continueParsing)
            {
                nextToken = token(is);

                if (!nextToken.isWord() || nextToken == token::END_SQR)
                {
                    continueParsing = false;
                }
            }
        }
        while (continueParsing);
    }
    else
    {
        // Read first five dimensions
        exponents_[dimensionSet::MASS] = nextToken.number();
        for (int d=1; d < dimensionSet::CURRENT; ++d)
        {
            is >> exponents_[d];
        }

        // Read next token
        token nextToken(is);

        // If next token is another number
        // read last two dimensions
        // and then read another token for the end of the dimensionSet
        if (nextToken.isNumber())
        {
            exponents_[dimensionSet::CURRENT] = nextToken.number();
            is >> nextToken;
            exponents_[dimensionSet::LUMINOUS_INTENSITY] = nextToken.number();
            is >> nextToken;
        }
        else
        {
            exponents_[dimensionSet::CURRENT] = 0;
            exponents_[dimensionSet::LUMINOUS_INTENSITY] = 0;
        }

        // Check end of dimensionSet
        if (nextToken != token::END_SQR)
        {
            FatalIOErrorInFunction(is)
                << "Expected a '" << token::END_SQR << "' in dimensionSet\n"
                << "in stream " << is.info() << nl
                << exit(FatalIOError);
        }
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::dimensionSet::write
(
    Ostream& os,
    scalar& multiplier,
    const dimensionSets& writeUnits
) const
{
    multiplier = 1.0;

    os << token::BEGIN_SQR;

    if (writeUnits.good() && os.format() == IOstreamOption::ASCII)
    {
        scalarField exponents(dimensionSet::nDimensions);
        for (int d=0; d < dimensionSet::nDimensions; ++d)
        {
            exponents[d] = exponents_[d];
        }
        writeUnits.coefficients(exponents);

        bool hasPrinted = false;

        // Set precision to lots
        std::streamsize oldPrecision = os.precision
        (
            std::numeric_limits<scalar>::digits10
        );

        forAll(exponents, i)
        {
            if (mag(exponents[i]) > smallExponent)
            {
                const dimensionedScalar& ds = writeUnits.units()[i];

                if (hasPrinted)
                {
                    os  << token::SPACE;
                }
                hasPrinted = true;
                os  << ds.name();
                if (mag(exponents[i]-1) > smallExponent)
                {
                    os  << '^' << exponents[i];

                    multiplier *= Foam::pow(ds.value(), exponents[i]);
                }
                else
                {
                    multiplier *= ds.value();
                }
            }
        }

        // Reset precision
        os.precision(oldPrecision);
    }
    else
    {
        for (int d=0; d < dimensionSet::nDimensions; ++d)
        {
            if (d) os << token::SPACE;
            os << exponents_[d];
        }
    }

    os  << token::END_SQR;

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Ostream& Foam::dimensionSet::write
(
    Ostream& os,
    scalar& multiplier
) const
{
    return write(os, multiplier, writeUnitSet());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, dimensionSet& ds)
{
    scalar mult(1.0);
    ds.read(is, mult);

    if (mag(mult-1.0) > dimensionSet::smallExponent)
    {
        FatalIOErrorInFunction(is)
            << "Cannot use scaled units in dimensionSet"
            << exit(FatalIOError);
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const dimensionSet& ds)
{
    scalar mult(1.0);
    ds.write(os, mult);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
