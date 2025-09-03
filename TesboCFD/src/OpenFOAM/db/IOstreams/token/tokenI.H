/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include <utility>  // For std::move, std::swap

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::token Foam::token::boolean(bool on) noexcept
{
    token tok;
    tok.type_ = tokenType::BOOL;
    tok.data_.labelVal = on;

    return tok;
}


inline Foam::token Foam::token::flag(int bitmask) noexcept
{
    token tok;
    tok.type_ = tokenType::FLAG;
    tok.data_.flagVal = bitmask;

    return tok;
}


inline bool Foam::token::is_wordToken(tokenType tokType) noexcept
{
    return
    (
        tokType == tokenType::WORD
     || tokType == tokenType::DIRECTIVE
    );
}


inline bool Foam::token::is_stringToken(tokenType tokType) noexcept
{
    return
    (
        tokType == tokenType::STRING
     || tokType == tokenType::EXPRESSION
     || tokType == tokenType::VARIABLE
     || tokType == tokenType::VERBATIM
     || tokType == tokenType::CHAR_DATA
    );
}


inline bool Foam::token::isseparator(int c) noexcept
{
    // NOTE: keep synchronized with ISstream::read(token&)

    switch (c)
    {
        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::PLUS :
        // Excluded token::MINUS since it could start a number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            return true;
        }

        default:
            break;
    }

    return false;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::token::setUndefined() noexcept
{
    type_ = tokenType::UNDEFINED;
    data_.int64Val = 0; // bit-wise zero for union content
    // leave lineNumber untouched - may still be needed
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::token::token() noexcept
:
    data_(), // bit-wise zero for union content
    type_(tokenType::UNDEFINED),
    line_(0)
{}


inline Foam::token::token(const token& tok)
:
    data_(tok.data_), // bit-wise copy of union content
    type_(tok.type_),
    line_(tok.line_)
{
    // Fundamental: values already handled by bit-wise copy
    // Pointer: duplicate content or increase refCount

    switch (type_)
    {
        // token::isWord()
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            data_.wordPtr = new word(*tok.data_.wordPtr);
            break;
        }

        // token::isString()
        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
        case tokenType::CHAR_DATA:
        {
            data_.stringPtr = new string(*tok.data_.stringPtr);
            break;
        }

        case tokenType::COMPOUND:
        {
            // Identical pointers, but increase the refCount
            data_.compoundPtr = tok.data_.compoundPtr;
            data_.compoundPtr->refCount::operator++();
            break;
        }

        default:
            break;
    }
}


inline Foam::token::token(token&& tok) noexcept
:
    data_(tok.data_), // bit-wise copy of union content
    type_(tok.type_),
    line_(tok.line_)
{
    tok.setUndefined(); // zero the union content without any checking
    tok.line_ = 0;
}


inline Foam::token::token(punctuationToken p, label lineNum) noexcept
:
    data_(),
    type_(tokenType::PUNCTUATION),
    line_(lineNum)
{
    data_.punctuationVal = p;
}


inline Foam::token::token(const label val, label lineNum) noexcept
:
    data_(),
    type_(tokenType::LABEL),
    line_(lineNum)
{
    data_.labelVal = val;
}


inline Foam::token::token(const float val, label lineNum) noexcept
:
    data_(),
    type_(tokenType::FLOAT),
    line_(lineNum)
{
    data_.floatVal = val;
}


inline Foam::token::token(const double val, label lineNum) noexcept
:
    data_(),
    type_(tokenType::DOUBLE),
    line_(lineNum)
{
    data_.doubleVal = val;
}


inline Foam::token::token(const word& w, label lineNum)
:
    data_(),
    type_(tokenType::WORD),
    line_(lineNum)
{
    data_.wordPtr = new word(w);
}


inline Foam::token::token(const string& str, label lineNum)
:
    data_(),
    type_(tokenType::STRING),
    line_(lineNum)
{
    data_.stringPtr = new string(str);
}


inline Foam::token::token(word&& w, label lineNum)
:
    data_(),
    type_(tokenType::WORD),
    line_(lineNum)
{
    data_.wordPtr = new word(std::move(w));
}


inline Foam::token::token(string&& str, label lineNum)
:
    data_(),
    type_(tokenType::STRING),
    line_(lineNum)
{
    data_.stringPtr = new string(std::move(str));
}


inline Foam::token::token
(
    tokenType tokType,
    const std::string& str,
    label lineNum
)
:
    data_(),
    type_(tokenType::STRING),
    line_(lineNum)
{
    if (is_wordToken(tokType))
    {
        type_ = tokType;
        data_.wordPtr = new word(str, false);  // no stripping
    }
    else
    {
        if (is_stringToken(tokType)) type_ = tokType;
        data_.stringPtr = new string(str);  // never strips
    }
}


inline Foam::token::token
(
    tokenType tokType,
    std::string&& str,
    label lineNum
)
:
    data_(),
    type_(tokenType::STRING),
    line_(lineNum)
{
    if (is_wordToken(tokType))
    {
        type_ = tokType;
        data_.wordPtr = new word(std::move(str), false);  // no stripping
    }
    else
    {
        if (is_stringToken(tokType)) type_ = tokType;
        data_.stringPtr = new string(std::move(str));  // never strips
    }
}


inline Foam::token::token(token::compound* ptr, label lineNum)
:
    data_(),
    type_(tokenType::COMPOUND),
    line_(lineNum)
{
    data_.compoundPtr = ptr;
    if (!data_.compoundPtr)
    {
        // Could handle nullptr as Fatal, but this is simpler
        setUndefined();
    }
}


inline Foam::token::token(autoPtr<token::compound>&& ptr, label lineNum)
:
    token(ptr.release(), lineNum)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::token::~token()
{
    reset();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::token::reset()
{
    switch (type_)
    {
        // token::isWord()
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            delete data_.wordPtr;
            break;
        }

        // token::isString()
        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
        case tokenType::CHAR_DATA:
        {
            delete data_.stringPtr;
            break;
        }

        case tokenType::COMPOUND:
        {
            if (data_.compoundPtr->refCount::unique())
            {
                delete data_.compoundPtr;
            }
            else
            {
                data_.compoundPtr->refCount::operator--();
            }
            break;
        }

        default:
            break;
    }

    setUndefined();
}


inline void Foam::token::setBad()
{
    reset();
    type_ = tokenType::ERROR;
}


inline void Foam::token::swap(token& tok) noexcept
{
    if (this == &tok)
    {
        return;  // Self-swap is a no-op
    }

    std::swap(data_, tok.data_);
    std::swap(type_, tok.type_);
    std::swap(line_, tok.line_);
}


inline Foam::word Foam::token::name() const
{
    return token::name(type_);
}


inline Foam::token::tokenType Foam::token::type() const noexcept
{
    return type_;
}


inline bool Foam::token::setType(token::tokenType tokType) noexcept
{
    if (type_ == tokType)
    {
        // No change required
        return true;
    }

    switch (tokType)
    {
        case tokenType::BOOL:
        case tokenType::LABEL:
        {
            switch (type_)
            {
                case tokenType::BOOL:
                case tokenType::LABEL:
                    type_ = tokType;
                    return true;
                    break;

                default:
                    break;
            }
        }
        break;

        // token::isWord()
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            switch (type_)
            {
                // token::isWord()
                case tokenType::WORD:
                case tokenType::DIRECTIVE:
                    type_ = tokType;
                    return true;
                    break;

                default:
                    break;
            }
        }
        break;

        // token::isString()
        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
        case tokenType::CHAR_DATA:
        {
            switch (type_)
            {
                // token::isWord()
                // could also go from WORD to STRING etc - to be decided

                // token::isString()
                case tokenType::STRING:
                case tokenType::EXPRESSION:
                case tokenType::VARIABLE:
                case tokenType::VERBATIM:
                case tokenType::CHAR_DATA:
                    type_ = tokType;
                    return true;
                    break;

                default:
                    break;
            }
        }
        break;

        default:
            break;
    }

    return false;
}


inline Foam::label Foam::token::lineNumber() const noexcept
{
    return line_;
}


inline Foam::label Foam::token::lineNumber(const label lineNum) noexcept
{
    label old(line_);
    line_ = lineNum;
    return old;
}


inline bool Foam::token::good() const noexcept
{
    return (type_ != tokenType::UNDEFINED && type_ != tokenType::ERROR);
}


inline bool Foam::token::undefined() const noexcept
{
    return (type_ == tokenType::UNDEFINED);
}


inline bool Foam::token::error() const noexcept
{
    return (type_ == tokenType::ERROR);
}


inline bool Foam::token::isBool() const noexcept
{
    return (type_ == tokenType::BOOL);
}


inline bool Foam::token::boolToken() const
{
    if (type_ == tokenType::BOOL || type_ == tokenType::LABEL)
    {
        return data_.labelVal;
    }

    parseError("bool");
    return false;
}


inline bool Foam::token::isFlag() const noexcept
{
    return (type_ == tokenType::FLAG);
}


inline int Foam::token::flagToken() const
{
    if (type_ == tokenType::FLAG)
    {
        return data_.flagVal;
    }

    parseError("flag bitmask");
    return flagType::NO_FLAG;
}


inline bool Foam::token::isPunctuation() const noexcept
{
    return (type_ == tokenType::PUNCTUATION);
}


inline bool Foam::token::isPunctuation(const punctuationToken p) const noexcept
{
    return
    (
        type_ == tokenType::PUNCTUATION
     && data_.punctuationVal == p
    );
}


inline bool Foam::token::isSeparator() const noexcept
{
    return
    (
        type_ == tokenType::PUNCTUATION
     && isseparator(data_.punctuationVal)
    );
}


inline Foam::token::punctuationToken Foam::token::pToken() const
{
    if (type_ == tokenType::PUNCTUATION)
    {
        return data_.punctuationVal;
    }

    parseError("punctuation character");
    return punctuationToken::NULL_TOKEN;
}


inline bool Foam::token::isLabel() const noexcept
{
    return
    (
        type_ == tokenType::LABEL
        // FUTURE?
        // || type_ == tokenType::INT32
        // || type_ == tokenType::INT64
    );
}


inline bool Foam::token::isLabel(const label value) const noexcept
{
    // FUTURE?
    // return
    // (
    //     type_ == tokenType::LABEL
    //   ? value == data_.labelVal
    //   : type_ == tokenType::INT32
    //   ? value == data_.int32Val
    //   : type_ == tokenType::INT64
    //   ? value == data_.int64Val
    //   : false
    // );

    return
    (
        type_ == tokenType::LABEL
     && value == data_.labelVal
    );
}


inline Foam::label Foam::token::labelToken() const
{
    if (type_ == tokenType::LABEL)
    {
        return data_.labelVal;
    }

    parseError("label");
    return 0;
}


inline bool Foam::token::isFloat() const noexcept
{
    return (type_ == tokenType::FLOAT);
}


inline float Foam::token::floatToken() const
{
    if (type_ == tokenType::FLOAT)
    {
        return data_.floatVal;
    }

    parseError("float");
    return 0;
}


inline bool Foam::token::isDouble() const noexcept
{
    return (type_ == tokenType::DOUBLE);
}


inline double Foam::token::doubleToken() const
{
    if (type_ == tokenType::DOUBLE)
    {
        return data_.doubleVal;
    }

    parseError("double");
    return 0;
}


inline bool Foam::token::isScalar() const noexcept
{
    return
    (
        type_ == tokenType::FLOAT
     || type_ == tokenType::DOUBLE
    );
}


inline Foam::scalar Foam::token::scalarToken() const
{
    if (type_ == tokenType::FLOAT)
    {
        return data_.floatVal;
    }
    else if (type_ == tokenType::DOUBLE)
    {
        return data_.doubleVal;
    }

    parseError("scalar");
    return 0;
}


inline bool Foam::token::isNumber() const noexcept
{
    return (isLabel() || isScalar());
}


inline Foam::scalar Foam::token::number() const
{
    if (isLabel())
    {
        return labelToken();
    }
    if (isScalar())
    {
        return scalarToken();
    }

    parseError("number (label or scalar)");
    return 0;
}


inline bool Foam::token::isWord() const noexcept
{
    return is_wordToken(type_);
}


inline bool Foam::token::isWord(const std::string& s) const
{
    return (isWord() && s == *data_.wordPtr);
}


inline bool Foam::token::isDirective() const noexcept
{
    return (type_ == tokenType::DIRECTIVE);
}


inline const Foam::word& Foam::token::wordToken() const
{
    if (isWord())
    {
        return *data_.wordPtr;
    }

    parseError("word");
    return word::null;
}


inline bool Foam::token::isQuotedString() const noexcept
{
    return (type_ == tokenType::STRING);
}


inline bool Foam::token::isString() const noexcept
{
    return is_stringToken(type_);
}


inline bool Foam::token::isExpression() const noexcept
{
    return (type_ == tokenType::EXPRESSION);
}


inline bool Foam::token::isVariable() const noexcept
{
    return (type_ == tokenType::VARIABLE);
}


inline bool Foam::token::isVerbatim() const noexcept
{
    return (type_ == tokenType::VERBATIM);
}


inline bool Foam::token::isCharData() const noexcept
{
    return (type_ == tokenType::CHAR_DATA);
}


inline bool Foam::token::isStringType() const noexcept
{
    return (isWord() || isString());
}


inline const Foam::string& Foam::token::stringToken() const
{
    if (isString())
    {
        return *data_.stringPtr;
    }
    else if (isWord())
    {
        // Foam::word derives from Foam::string, no need to cast.
        return *data_.wordPtr;
    }

    parseError("string");
    return string::null;
}


// Could also have reference to stringToken() - eg refStringToken() ?
// but not really sure where this would be useful...

inline bool Foam::token::isCompound() const noexcept
{
    return (type_ == tokenType::COMPOUND);
}


inline bool Foam::token::isCompound(const word& compoundType) const
{
    return
    (
        type_ == tokenType::COMPOUND
     && data_.compoundPtr->type() == compoundType
    );
}


template<class Type>
inline const Type* Foam::token::isCompound() const
{
    return
    (
        type_ == tokenType::COMPOUND
      ? dynamic_cast<const Type*>(data_.compoundPtr)
      : nullptr
    );
}


inline const Foam::token::compound& Foam::token::compoundToken() const
{
    if (type_ != tokenType::COMPOUND)
    {
        parseError("compound");
    }
    return *data_.compoundPtr;
}


inline Foam::token::compound& Foam::token::refCompoundToken()
{
    if (type_ != tokenType::COMPOUND)
    {
        parseError("compound");
    }
    return *data_.compoundPtr;
}


template<class Type>
inline Type& Foam::token::refCompoundToken()
{
    if (type_ != tokenType::COMPOUND)
    {
        parseError("compound");
    }
    return static_cast<Type&>
    (
        dynamicCast<token::Compound<Type>>
        (
            *data_.compoundPtr
        )
    );
}


inline Foam::token::compound&
Foam::token::transferCompoundToken(const Istream& is)
{
    return transferCompoundToken(&is);
}


template<class Type>
inline Type& Foam::token::transferCompoundToken(const Istream& is)
{
    return static_cast<Type&>
    (
        dynamicCast<token::Compound<Type>>
        (
            transferCompoundToken(&is)
        )
    );
}


template<class T>
Foam::token::tokenType Foam::token::Compound<T>::typeCode() const
{
    // Does not cover all possibilities perfectly, but should handle
    // most of the common ones (bool, label, scalar, vector lists).
    // Something like List<edge> will not be quite correct if we rely
    // on nComponents

    typedef typename T::value_type valueType;

    if (std::is_same<valueType, bool>::value)
    {
        // List<bool>
        return token::tokenType::BOOL;
    }
    else if (is_contiguous_label<valueType>::value)
    {
        // List<label>, List<labelVector> etc
        return token::tokenType::LABEL;

        // FUTURE?
        // return
        // (
        //     sizeof(int32_t) == sizeof(Foam::label)
        //   ? token::tokenType::INT32
        //   : token::tokenType::INT64
        // );

    }
    else if (is_contiguous_scalar<valueType>::value)
    {
        // List<scalar>, List<vector>, List<tensor> etc
        return
        (
            sizeof(float) == sizeof(Foam::scalar)
          ? token::tokenType::FLOAT
          : token::tokenType::DOUBLE
        );
    }
    else if (std::is_same<valueType, char>::value)
    {
        // List<char>
        return token::tokenType::PUNCTUATION;
    }
    else
    {
        // Do not handle List<word> or List<string> at the moment
        // since filling non-contiguous data is probably not desirable
        return token::tokenType::UNDEFINED;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::token::operator=(const token& tok)
{
    if (this == &tok)
    {
        return;  // Self-assignment is a no-op
    }

    reset();

    type_ = tok.type_;
    data_ = tok.data_;  // bit-wise copy of union content
    line_ = tok.line_;

    // Fundamental: values already handled by bit-wise copy
    // Pointer: duplicate content or increase refCount

    switch (type_)
    {
        // token::isWord()
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            data_.wordPtr = new word(*tok.data_.wordPtr);
        }
        break;

        // token::isString()
        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
        case tokenType::CHAR_DATA:
        {
            data_.stringPtr = new string(*tok.data_.stringPtr);
        }
        break;

        case tokenType::COMPOUND:
        {
            // Identical pointers, but increase the refCount
            data_.compoundPtr = tok.data_.compoundPtr;
            data_.compoundPtr->refCount::operator++();
        }
        break;

        default:
            break;
    }
}


inline void Foam::token::operator=(token&& tok)
{
    if (this == &tok)
    {
        return;  // Self-assignment is a no-op
    }

    reset();
    line_ = 0;
    swap(tok);
}


inline void Foam::token::operator=(const punctuationToken p)
{
    reset();
    type_ = tokenType::PUNCTUATION;
    data_.punctuationVal = p;
}


inline void Foam::token::operator=(const label val)
{
    reset();
    type_ = tokenType::LABEL;
    data_.labelVal = val;
}


inline void Foam::token::operator=(const float val)
{
    reset();
    type_ = tokenType::FLOAT;
    data_.floatVal = val;
}


inline void Foam::token::operator=(const double val)
{
    reset();
    type_ = tokenType::DOUBLE;
    data_.doubleVal = val;
}


inline void Foam::token::operator=(const word& w)
{
    reset();
    type_ = tokenType::WORD;
    data_.wordPtr = new word(w);
}


inline void Foam::token::operator=(const string& str)
{
    reset();
    type_ = tokenType::STRING;
    data_.stringPtr = new string(str);
}


inline void Foam::token::operator=(word&& w)
{
    reset();
    type_ = tokenType::WORD;
    data_.wordPtr = new word(std::move(w));
}


inline void Foam::token::operator=(string&& s)
{
    reset();
    type_ = tokenType::STRING;
    data_.stringPtr = new string(std::move(s));
}


inline void Foam::token::operator=(Foam::token::compound* ptr)
{
    reset();
    type_ = tokenType::COMPOUND;
    data_.compoundPtr = ptr;
}


inline void Foam::token::operator=(autoPtr<token::compound>&& ptr)
{
    reset();
    type_ = tokenType::COMPOUND;
    data_.compoundPtr = ptr.release();
}


inline bool Foam::token::operator==(const token& tok) const
{
    if (type_ != tok.type_)
    {
        return false;
    }

    switch (type_)
    {
        case tokenType::UNDEFINED:
            return true;

        case tokenType::BOOL:
            return data_.labelVal == tok.data_.labelVal;

        case tokenType::FLAG:
            return data_.flagVal == tok.data_.flagVal;

        case tokenType::PUNCTUATION:
            return data_.punctuationVal == tok.data_.punctuationVal;

        case tokenType::LABEL:
            return data_.labelVal == tok.data_.labelVal;

        case tokenType::FLOAT:
            return equal(data_.floatVal, tok.data_.floatVal);

        case tokenType::DOUBLE:
            return equal(data_.doubleVal, tok.data_.doubleVal);

        // token::isWord()
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
            return *data_.wordPtr == *tok.data_.wordPtr;

        // token::isString()
        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
        case tokenType::CHAR_DATA:
            return *data_.stringPtr == *tok.data_.stringPtr;

        case tokenType::COMPOUND:
            return data_.compoundPtr == tok.data_.compoundPtr;

        case tokenType::ERROR:
            return true;
    }

    return false;
}


inline bool Foam::token::operator==(const punctuationToken p) const noexcept
{
    return isPunctuation(p);
}


inline bool Foam::token::operator==(const std::string& s) const
{
    return
    (
        isWord()
      ? s == *data_.wordPtr
      : isString() && s == *data_.stringPtr
    );
}


inline bool Foam::token::operator==(const label val) const noexcept
{
    return isLabel(val);
}


inline bool Foam::token::operator==(const float val) const noexcept
{
    return
    (
        type_ == tokenType::FLOAT
     && equal(data_.floatVal, val)
    );
}


inline bool Foam::token::operator==(const double val) const noexcept
{
    return
    (
        type_ == tokenType::DOUBLE
     && equal(data_.doubleVal, val)
    );
}


inline bool Foam::token::operator!=(const token& tok) const
{
    return !operator==(tok);
}


inline bool Foam::token::operator!=(const punctuationToken p) const noexcept
{
    return !isPunctuation(p);
}


inline bool Foam::token::operator!=(const label val) const noexcept
{
    return !operator==(val);
}


inline bool Foam::token::operator!=(const float val) const noexcept
{
    return !operator==(val);
}


inline bool Foam::token::operator!=(const double val) const noexcept
{
    return !operator==(val);
}


inline bool Foam::token::operator!=(const std::string& s) const
{
    return !operator==(s);
}


// ************************************************************************* //
