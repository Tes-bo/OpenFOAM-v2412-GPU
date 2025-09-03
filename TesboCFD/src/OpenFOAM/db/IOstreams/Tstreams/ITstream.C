/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "error.H"
#include "ITstream.H"
#include "SpanStream.H"
#include <algorithm>
#include <memory>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

static std::unique_ptr<Foam::ITstream> emptyStreamPtr_;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Convert input sequence into a list of tokens.
// Return the number of tokens in the resulting list.
static label parseStream(ISstream& is, tokenList& tokens)
{
    tokens.clear();

    label count = 0;
    token tok;
    while (!is.read(tok).bad() && tok.good())
    {
        if (count >= tokens.size())
        {
            // Increase capacity (doubling) with min-size [64]
            tokens.resize(Foam::max(label(64), label(2*tokens.size())));
        }

        tokens[count] = std::move(tok);
        ++count;
    }

    tokens.resize(count);

    return count;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::ITstream& Foam::ITstream::empty_stream()
{
    if (emptyStreamPtr_)
    {
        emptyStreamPtr_->ITstream::clear();  // Ensure it really is empty
        emptyStreamPtr_->ITstream::seek(0);  // rewind() bypassing virtual
    }
    else
    {
        emptyStreamPtr_.reset(new ITstream(Foam::zero{}, "empty-stream"));
    }

    // Set stream as bad to indicate that this is an invald stream
    emptyStreamPtr_->setBad();

    return *emptyStreamPtr_;
}


Foam::tokenList Foam::ITstream::parse_chars
(
    const char* s,
    size_t nbytes,
    IOstreamOption streamOpt
)
{
    ISpanStream is(s, nbytes, streamOpt);

    tokenList tokens;
    parseStream(is, tokens);
    return tokens;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ITstream::reset(const char* input, size_t nbytes)
{
    ISpanStream is(input, nbytes, static_cast<IOstreamOption>(*this));

    parseStream(is, static_cast<tokenList&>(*this));
    ITstream::seek(0);  // rewind() bypassing virtual
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ITstream::reserveCapacity(const label newCapacity)
{
    // Reserve - leave excess capacity for further appends

    label len = tokenList::size();

    if (len < newCapacity)
    {
        // Min-size (16) when starting from zero
        if (!len) len = 8;

        // Increase capacity. Strict doubling
        do
        {
            len *= 2;
        }
        while (len < newCapacity);

        tokenList::resize(len);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ITstream::ITstream(const ITstream& is)
:
    Istream(static_cast<IOstreamOption>(is)),
    tokenList(is),
    name_(is.name_),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream(ITstream&& is)
:
    Istream(static_cast<IOstreamOption>(is)),
    tokenList(std::move(static_cast<tokenList&>(is))),
    name_(std::move(is.name_)),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    IOstreamOption streamOpt,
    const string& name
)
:
    Istream(IOstreamOption(streamOpt.format(), streamOpt.version())),
    tokenList(),
    name_(name),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    const Foam::zero,
    const string& name,
    IOstreamOption streamOpt
)
:
    ITstream(streamOpt, name)
{}


Foam::ITstream::ITstream
(
    const UList<token>& tokens,
    IOstreamOption streamOpt,
    const string& name
)
:
    Istream(IOstreamOption(streamOpt.format(), streamOpt.version())),
    tokenList(tokens),
    name_(name),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    List<token>&& tokens,
    IOstreamOption streamOpt,
    const string& name
)
:
    Istream(IOstreamOption(streamOpt.format(), streamOpt.version())),
    tokenList(std::move(tokens)),
    name_(name),
    tokenIndex_(0)
{
    setOpened();
    setGood();
}


Foam::ITstream::ITstream
(
    const UList<char>& input,
    IOstreamOption streamOpt,
    const string& name
)
:
    ITstream(streamOpt, name)
{
    reset(input.cdata(), input.size_bytes());
}


Foam::ITstream::ITstream
(
    const std::string& input,
    IOstreamOption streamOpt,
    const string& name
)
:
    ITstream(streamOpt, name)
{
    reset(input.data(), input.size());
}


Foam::ITstream::ITstream
(
    const char* input,
    IOstreamOption streamOpt,
    const string& name
)
:
    ITstream(streamOpt, name)
{
    reset(input, strlen(input));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ITstream::print(Ostream& os) const
{
    os  << "ITstream : " << name_.c_str() << ", line ";

    if (tokenList::empty())
    {
        os  << lineNumber();
    }
    else
    {
        const tokenList& toks = *this;

        os  << toks.front().lineNumber();

        if (toks.front().lineNumber() < toks.back().lineNumber())
        {
            os  << '-' << toks.back().lineNumber();
        }
    }
    os  << ", ";

    IOstream::print(os);
}


std::string Foam::ITstream::toString() const
{
    if (tokenList::empty())
    {
        return std::string();
    }
    else if (tokenList::size() == 1 && tokenList::front().isStringType())
    {
        // Already a string-type (WORD, STRING, ...). Just copy.
        return tokenList::front().stringToken();
    }

    // Stringify
    OCharStream buf;
    buf.precision(16);      // Some reasonably high precision

    auto iter = tokenList::cbegin();
    const auto last = tokenList::cend();

    // Note: could also just use the buffer token-wise

    // Contents - space separated
    if (iter != last)
    {
        buf << *iter;

        for (++iter; (iter != last); (void)++iter)
        {
            buf << token::SPACE << *iter;
        }
    }

    const auto view = buf.view();

    return std::string(view.data(), view.size());
}


const Foam::token& Foam::ITstream::peek() const noexcept
{
    // Use putback token if it exists
    if (Istream::hasPutback())
    {
        return Istream::peekBack();
    }

    return peekNoFail(tokenIndex_);
}


Foam::token& Foam::ITstream::currentToken()
{
    if (tokenIndex_ < 0 || tokenIndex_ >= tokenList::size())
    {
        FatalIOErrorInFunction(*this)
            << "Token index " << tokenIndex_ << " out of range [0,"
            << tokenList::size() << "]\n"
            << abort(FatalIOError);
    }

    return tokenList::operator[](tokenIndex_);
}


void Foam::ITstream::seek(label pos) noexcept
{
    lineNumber_ = 0;

    if (pos < 0 || pos >= tokenList::size())
    {
        // Seek end (-1) or seek is out of range
        tokenIndex_ = tokenList::size();

        if (!tokenList::empty())
        {
            // The closest reference lineNumber
            lineNumber_ = tokenList::cdata()[tokenIndex_-1].lineNumber();
        }

        setEof();
    }
    else
    {
        tokenIndex_ = pos;

        if (tokenIndex_ < tokenList::size())
        {
            lineNumber_ = tokenList::cdata()[tokenIndex_].lineNumber();
        }

        setOpened();
        setGood();
    }
}


bool Foam::ITstream::skip(label n) noexcept
{
    if (!n)
    {
        // No movement - just check the current range
        return (tokenIndex_ >= 0 && tokenIndex_ < tokenList::size());
    }

    tokenIndex_ += n;  // Move forward (+ve) or backwards (-ve)

    bool noError = true;

    if (tokenIndex_ < 0)
    {
        // Underflow range
        noError = false;
        tokenIndex_ = 0;
    }
    else if (tokenIndex_ >= tokenList::size())
    {
        // Overflow range
        noError = false;
        tokenIndex_ = tokenList::size();

        if (!tokenList::empty())
        {
            // The closest reference lineNumber
            lineNumber_ = tokenList::cdata()[tokenIndex_-1].lineNumber();
        }
    }

    // Update stream information
    if (tokenIndex_ < tokenList::size())
    {
        lineNumber_ = tokenList::cdata()[tokenIndex_].lineNumber();
        setOpened();
        setGood();
    }
    else
    {
        setEof();
    }

    return noError;
}


Foam::Istream& Foam::ITstream::read(token& tok)
{
    // Use putback token if it exists
    if (Istream::getBack(tok))
    {
        lineNumber_ = tok.lineNumber();
        return *this;
    }

    tokenList& toks = *this;
    const label nToks = toks.size();

    if (tokenIndex_ < nToks)
    {
        tok = toks[tokenIndex_++];
        lineNumber_ = tok.lineNumber();

        if (tokenIndex_ == nToks)
        {
            setEof();
        }
    }
    else
    {
        if (eof())
        {
            FatalIOErrorInFunction(*this)
                << "attempt to read beyond EOF"
                << exit(FatalIOError);
            setBad();
        }
        else
        {
            setEof();
        }

        tok.reset();

        if (nToks)
        {
            tok.lineNumber(toks.back().lineNumber());
        }
        else
        {
            tok.lineNumber(this->lineNumber());
        }
    }

    return *this;
}


Foam::labelRange Foam::ITstream::find
(
    const token::punctuationToken delimOpen,
    const token::punctuationToken delimClose,
    label pos
) const
{
    if (pos < 0)
    {
        pos = tokenIndex_;
    }

    labelRange slice;

    for (label depth = 0; pos < tokenList::size(); ++pos)
    {
        const token& tok = tokenList::operator[](pos);

        if (tok.isPunctuation())
        {
            if (tok.isPunctuation(delimOpen))
            {
                if (!depth)
                {
                    // Initial open delimiter
                    slice.start() = pos;
                }

                ++depth;
            }
            else if (tok.isPunctuation(delimClose))
            {
                --depth;

                if (depth < 0)
                {
                    // A closing delimiter without an open!
                    // Raise error?
                    break;
                }
                if (!depth)
                {
                    // The end - include delimiter into the count
                    slice.size() = (pos - slice.start()) + 1;
                    break;
                }
            }
        }
    }

    return slice;
}


Foam::ITstream Foam::ITstream::extract(const labelRange& range)
{
    ITstream result
    (
        static_cast<IOstreamOption>(*this),
        this->name()
    );
    result.setLabelByteSize(this->labelByteSize());
    result.setScalarByteSize(this->scalarByteSize());

    // Validate the slice range of list
    const labelRange slice(range.subset0(tokenList::size()));

    if (!slice.good())
    {
        // No-op
        return result;
    }

    auto first = tokenList::begin(slice.begin_value());
    auto last = tokenList::begin(slice.end_value());

    result.resize(label(last - first));

    // Move tokens into result list
    std::move(first, last, result.begin());
    result.seek(0);  // rewind() bypassing virtual


    (void) remove(slice);  // Adjust the original list

    return result;
}


Foam::label Foam::ITstream::remove(const labelRange& range)
{
    // Validate the slice range of list
    const labelRange slice(range.subset0(tokenList::size()));

    if (!slice.good())
    {
        // No-op
        return 0;
    }

    if (slice.end_value() >= tokenList::size())
    {
        // Remove entire tail
        tokenList::resize(slice.begin_value());
    }
    else
    {
        // Attempt to adjust the current token index to something sensible...
        if (slice.contains(tokenIndex_))
        {
            // Within the removed slice - reposition tokenIndex before it
            seek(slice.begin_value());
            skip(-1);
        }
        else if (tokenIndex_ >= slice.end_value())
        {
            // After the removed slice - reposition tokenIndex relatively
            skip(-slice.size());
        }

        // Move tokens down in the list
        std::move
        (
            tokenList::begin(slice.end_value()),
            tokenList::end(),
            tokenList::begin(slice.begin_value())
        );

        // Truncate
        tokenList::resize(tokenList::size() - slice.size());
    }

    if (tokenIndex_ >= tokenList::size())
    {
        tokenIndex_ = tokenList::size();
        setEof();
    }
    else if (tokenIndex_ >= 0 && tokenIndex_ < tokenList::size())
    {
        lineNumber_ = tokenList::operator[](tokenIndex_).lineNumber();
    }

    return slice.size();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::ITstream::read(char&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(word&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(string&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(label&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(float&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(double&)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::readRaw(char*, std::streamsize)
{
    NotImplemented;
    return *this;
}


Foam::Istream& Foam::ITstream::read(char*, std::streamsize)
{
    NotImplemented;
    return *this;
}


void Foam::ITstream::add_tokens(const token& tok)
{
    reserveCapacity(tokenIndex_ + 1);

    tokenList::operator[](tokenIndex_) = tok;
    ++tokenIndex_;
}


void Foam::ITstream::add_tokens(token&& tok)
{
    reserveCapacity(tokenIndex_ + 1);

    tokenList::operator[](tokenIndex_) = std::move(tok);
    ++tokenIndex_;
}


void Foam::ITstream::add_tokens(const UList<token>& toks)
{
    const label len = toks.size();
    reserveCapacity(tokenIndex_ + len);

    std::copy_n(toks.begin(), len, tokenList::begin(tokenIndex_));
    tokenIndex_ += len;
}


void Foam::ITstream::add_tokens(List<token>&& toks)
{
    const label len = toks.size();
    reserveCapacity(tokenIndex_ + len);

    std::move(toks.begin(), toks.end(), tokenList::begin(tokenIndex_));
    tokenIndex_ += len;
    toks.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::ITstream::operator=(const ITstream& is)
{
    // Self-assignment is a no-op
    if (this != &is)
    {
        Istream::operator=(is);
        tokenList::operator=(is);
        name_ = is.name_;
        ITstream::seek(0);  // rewind() bypassing virtual
    }
}


void Foam::ITstream::operator=(const UList<token>& toks)
{
    tokenList::operator=(toks);
    ITstream::seek(0);  // rewind() bypassing virtual
}


void Foam::ITstream::operator=(List<token>&& toks)
{
    tokenList::operator=(std::move(toks));
    ITstream::seek(0);  // rewind() bypassing virtual
}


// ************************************************************************* //
