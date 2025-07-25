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

Class
    Foam::token

Description
    A token holds an item read from Istream.

SourceFiles
    tokenI.H
    token.C
    tokenIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_token_H
#define Foam_token_H

#include "label.H"
#include "uLabel.H"
#include "scalar.H"
#include "word.H"
#include "contiguous.H"
#include "refCount.H"
#include "InfoProxy.H"
#include "typeInfo.H"

#define NoHashTableC
#include "runTimeSelectionTables.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class token;
Ostream& operator<<(Ostream& os, const token& tok);

/*---------------------------------------------------------------------------*\
                            Class token Declaration
\*---------------------------------------------------------------------------*/

class token
{
public:

    //- Enumeration defining the types of token.
    //  Since the enumeration is used to tag content in Pstream, it is of
    //  type \c char and shall have values that do not overlap with regular
    //  punctuation characters.
    enum tokenType : char
    {
        UNDEFINED = '\0',     //!< An undefined token-type
        ERROR     = '\x80',   //!< Token error encountered

        // Fundamental types
        FLAG,             //!< stream flag (1-byte bitmask)
        PUNCTUATION,      //!< single character punctuation
        BOOL,             //!< boolean type
        LABEL,            //!< label (integer) type
        FLOAT,            //!< float (single-precision) type
        DOUBLE,           //!< double (double-precision) type

        // Pointer types
        WORD,             //!< Foam::word
        STRING,           //!< Foam::string (usually double-quoted)
        COMPOUND,         //!< Compound type such as \c List\<label\> etc.

        DIRECTIVE,        //!< Word-variant: dictionary \c \#directive
                          //!< stored with sigil

        EXPRESSION,       //!< String-variant: math expression for evaluation
                          //!< stored with delimiters
        VARIABLE,         //!< String-variant: dictionary \c \$variable
                          //!< stored with sigil
        VERBATIM,         //!< String-variant: verbatim string content
                          //!< stored without delimiters
        CHAR_DATA,        //!< String-variant: plain character content

        // Aliases
        FLOAT_SCALAR = FLOAT,       //!< compatibility name for FLOAT
        DOUBLE_SCALAR = DOUBLE,     //!< compatibility name for DOUBLE
        VERBATIMSTRING = VERBATIM   //!< compatibility name for VERBATIM
    };


    //- Stream or output control flags (1-byte width)
    enum flagType
    {
        NO_FLAG = 0,   //!< No flags
        ASCII   = 1,   //!< ASCII-mode stream
        BINARY  = 2    //!< BINARY-mode stream
    };


    //- Standard punctuation tokens (a character)
    enum punctuationToken : char
    {
        NULL_TOKEN     = '\0',  //!< Nul character
        TAB            = '\t',  //!< Tab [isspace]
        NL             = '\n',  //!< Newline [isspace]
        SPACE          = ' ',   //!< Space [isspace]

        COLON          = ':',   //!< Colon [#isseparator]
        SEMICOLON      = ';',   //!< Semicolon [#isseparator]
        COMMA          = ',',   //!< Comma [#isseparator]
        HASH           = '#',   //!< Hash - directive or start verbatim string
        DOLLAR         = '$',   //!< Dollar - start variable or expression
        QUESTION       = '?',   //!< Question mark (eg, ternary)
        ATSYM          = '@',   //!< The 'at' symbol
        SQUOTE         = '\'',  //!< Single quote
        DQUOTE         = '"',   //!< Double quote

        ASSIGN         = '=',   //!< Assignment/equals [#isseparator]
        PLUS           = '+',   //!< Addition [#isseparator]
        MINUS          = '-',   //!< Subtract or start of negative number
        MULTIPLY       = '*',   //!< Multiply [#isseparator]
        DIVIDE         = '/',   //!< Divide [#isseparator]

        LPAREN         = '(',   //!< Left parenthesis [#isseparator]
        RPAREN         = ')',   //!< Right parenthesis [#isseparator]
        LSQUARE        = '[',   //!< Left square bracket [#isseparator]
        RSQUARE        = ']',   //!< Right square bracket [#isseparator]
        LBRACE         = '{',   //!< Left brace [#isseparator]
        RBRACE         = '}',   //!< Right brace [#isseparator]

        // With semantically meaning

        ADD            = PLUS,      //!< Addition [#isseparator]
        SUBTRACT       = MINUS,     //!< Subtract or start of negative number
        END_STATEMENT  = SEMICOLON, //!< End entry [#isseparator]
        BEGIN_LIST     = LPAREN,    //!< Begin list [#isseparator]
        END_LIST       = RPAREN,    //!< End list [#isseparator]
        BEGIN_SQR      = LSQUARE,   //!< Begin dimensions [#isseparator]
        END_SQR        = RSQUARE,   //!< End dimensions [#isseparator]
        BEGIN_BLOCK    = LBRACE,    //!< Begin block [#isseparator]
        END_BLOCK      = RBRACE,    //!< End block [#isseparator]

        BEGIN_STRING   = DQUOTE,    //!< Begin string with double quote
        END_STRING     = DQUOTE     //!< End string with double quote
    };


    //- Abstract base class for complex tokens
    class compound
    :
        public refCount
    {
    protected:

        //- The compound token state, internally used values only
        //  (0: initial, 1: moved, 2: pending read...)
        unsigned char state_;

    public:

        //- Declare type-name, virtual type (without debug switch)
        TypeNameNoDebug("compound");

        //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            compound,
            empty,
            (),
            ()
        );


        // Constructors

            //- Default construct
            constexpr compound() noexcept
            :
                state_(0)
            {}


        //- Destructor
        virtual ~compound() noexcept = default;


        // Selectors

            //- Default construct specified compound type
            static autoPtr<compound> New(const word& compoundType);

            //- Default construct specified compound type
            //- and read from stream
            static autoPtr<compound> New
            (
                const word& compoundType,
                Istream& is,
                const bool readContent = true
            );

            //- Attempt dynamic_cast to \c Type
            //- returns nullptr if cast is not possible
            template<class Type>
            const Type* isA() const
            {
                return dynamic_cast<const Type*>(this);
            }


        // Member Functions

            //- True if a known (registered) compound type
            static bool isCompound(const word& compoundType);

            //- Get compound transferred status
            bool moved() const noexcept { return (state_ & 1); }

            //- Set compound transferred status
            void moved(bool b) noexcept
            {
                if (b) { state_ |= 1; } else { state_ &= ~1; }
            }

            //- Get compound pending-read status
            bool pending() const noexcept { return (state_ & 2); }

            //- Set compound pending-read status
            void pending(bool b) noexcept
            {
                if (b) { state_ |= 2; } else { state_ &= ~2; }
            }

            //- The size of the underlying content
            virtual label size() const = 0;

            //- Change the size of the underlying container content
            virtual void resize(const label n) = 0;

            //- Fill with zero values, or with default constructed
            virtual void fill_zero() = 0;

            //- Read from stream into the underlying content
            virtual void read(Istream& is) = 0;

            //- Write the underlying content
            virtual void write(Ostream& os) const = 0;


        // Attributes and access methods

            //- The token type (if any) corresponding to the contained
            //- component type (LABEL, FLOAT, DOUBLE, etc).
            virtual tokenType typeCode() const = 0;

            //- The number of vector-space or other components
            //- of the underlying data content
            virtual direction nComponents() const = 0;

            //- The vector-space rank associated with the data content
            virtual direction rank() const = 0;

            //- Pointer to the underlying data as byte data
            virtual const char* cdata_bytes() const = 0;

            //- Pointer to the underlying data as byte data
            virtual char* data_bytes() = 0;

            //- Size of the (contiguous) byte data
            virtual std::streamsize size_bytes() const = 0;


        // Operators

            //- Output operator
            friend Ostream& operator<<(Ostream& os, const compound& ct);
    };


    //- A templated class for holding compound tokens.
    //- The underlying container is normally a List of values,
    //- it must have a \c value_type typedef as well as
    //- size(), resize(), cdata_bytes(), data_bytes(), size_bytes() methods
    template<class T>
    class Compound
    :
        public token::compound,
        public T
    {
        // Private Member Functions

            //- Fill with zero (contiguous types)
            template<class ValT = typename T::value_type>
            typename std::enable_if<pTraits_has_zero<ValT>::value, void>::type
            _m_fill_zero()
            {
                T::operator=(pTraits<ValT>::zero);
            }

            //- Fill with default value initialized (non-contiguous types)
            template<class ValT = typename T::value_type>
            typename std::enable_if<!pTraits_has_zero<ValT>::value, void>::type
            _m_fill_zero()
            {
                T::operator=(ValT());
            }

    public:

        //- The type of values held by the compound
        typedef typename T::value_type value_type;

        //- Declare type-name, virtual type (without debug switch)
        TypeNameNoDebug("Compound<T>");

        //- No copy construct
        Compound(const Compound<T>&) = delete;

        //- No copy assignment
        void operator=(const Compound<T>&) = delete;

        // Constructors

            //- Default construct
            Compound() = default;

            //- Copy construct from base content
            explicit Compound(const T& content)
            :
                T(content)
            {}

            //- Move construct from base content
            explicit Compound(T&& content)
            :
                T(std::move(content))
            {}

            //- Read construct from Istream
            explicit Compound(Istream& is)
            :
                T(is)
            {}


        // Selectors

            //- Construct autoPtr compound with forwarding arguments.
            //  The return type is compound, not Compound since that
            //  is what the token interface requires.
            template<class... Args>
            static autoPtr<compound> New(Args&&... args)
            {
                return autoPtr<compound>
                (
                    new Compound<T>(T(std::forward<Args>(args)...))
                );
            }


        // Member Functions

            //- The size of the underlying content
            virtual label size() const
            {
                return T::size();
            }

            //- Change the size of the underlying content
            virtual void resize(const label n)
            {
                T::resize(n);
            }

            //- Fill with zero value or with default value initialized
            virtual void fill_zero() { _m_fill_zero(); }

            //- Redirect read to underlying content
            virtual void read(Istream& is)
            {
                token::compound::state_ = 0;  // Clean state
                operator>>(is, static_cast<T&>(*this));
            }

            //- Redirect write to underlying content
            virtual void write(Ostream& os) const
            {
                operator<<(os, static_cast<const T&>(*this));
            }


        // Attributes and access methods

            //- The token type (if any) corresponding to the contained
            //- component type (LABEL, FLOAT, DOUBLE, etc).
            virtual tokenType typeCode() const;

            //- The number of vector-space or other components
            //- of the underlying data content
            virtual direction nComponents() const
            {
                return pTraits_nComponents<value_type>();
            }

            //- The vector-space rank associated with the data content
            virtual direction rank() const
            {
                return pTraits_rank<value_type>();
            }

            //- Pointer to the underlying data as byte data
            virtual const char* cdata_bytes() const
            {
                return T::cdata_bytes();
            }

            //- Pointer to the underlying data as byte data
            virtual char* data_bytes()
            {
                return T::data_bytes();
            }

            //- Size of the (contiguous) byte data
            virtual std::streamsize size_bytes() const
            {
                return T::size_bytes();
            }
    };


    //- An undefined token
    static const token undefinedToken;


private:

    //- A %union of token types
    union content
    {
        // Fundamental values. Largest first for any {} initialization.
        int64_t int64Val;
        int32_t int32Val;

        int flagVal;   // bitmask - stored as int, not enum
        punctuationToken punctuationVal;
        label labelVal;
        float floatVal;
        double doubleVal;

        // Pointers
        word* wordPtr;
        string* stringPtr;
        mutable compound* compoundPtr;
    };


    // Private Data

        //- The data content (as a union).
        //  For memory alignment this should appear as the first member.
        content data_;

        //- The token type
        tokenType type_;

        //- The file line number where the token was read from
        label line_;


    // Private Member Functions

        //- Set as UNDEFINED and zero the union content without any checking
        inline void setUndefined() noexcept;

        //- True if token type corresponds to WORD or WORD-variant.
        inline static bool is_wordToken(tokenType tokType) noexcept;

        //- True if token type corresponds to STRING or STRING-variant
        inline static bool is_stringToken(tokenType tokType) noexcept;

        // Parse error, expected 'expected', found ...
        void parseError(const char* expected) const;


public:

    // Static Data Members

        //- The type name is "token"
        static constexpr const char* const typeName = "token";


    // Constructors

        //- Default construct, initialized to an UNDEFINED token.
        inline constexpr token() noexcept;

        //- Copy construct
        inline token(const token& t);

        //- Move construct. The original token is left as UNDEFINED.
        inline token(token&& t) noexcept;

        //- Construct punctuation character token
        inline explicit token(punctuationToken p, label lineNum=0) noexcept;

        //- Construct label token
        //- \note Use boolean() static method for creating a \b bool token
        //- since \c bool and \c int are otherwise indistinguishable
        inline explicit token(const label val, label lineNum=0) noexcept;

        //- Construct float token
        inline explicit token(const float val, label lineNum=0) noexcept;

        //- Construct double token
        inline explicit token(const double val, label lineNum=0) noexcept;

        //- Copy construct word token
        inline explicit token(const word& w, label lineNum=0);

        //- Copy construct string token
        inline explicit token(const string& str, label lineNum=0);

        //- Move construct word token
        inline explicit token(word&& w, label lineNum=0);

        //- Move construct string token
        inline explicit token(string&& str, label lineNum=0);

        //- Copy construct word/string token with the specified variant.
        //  A invalid word/string variant type is silently treated as STRING.
        //  No character stripping
        inline explicit token(tokenType typ, const std::string&, label line=0);

        //- Move construct word/string token with the specified variant.
        //  A invalid word/string variant type is silently treated as STRING.
        //  No character stripping
        inline explicit token(tokenType typ, std::string&&, label line=0);

        //- Construct from a compound pointer, taking ownership
        inline explicit token(token::compound* ptr, label lineNum=0);

        //- Move construct from a compound pointer, taking ownership
        inline explicit token(autoPtr<token::compound>&& ptr, label lineNum=0);

        //- Construct from Istream
        explicit token(Istream& is);


    //- Destructor
    inline ~token();


    // Static Functions

        //- Create a bool token.
        inline static token boolean(bool on) noexcept;

        //- Create a token with stream flags, no sanity check
        //
        //  \param bitmask the flags to set
        inline static token flag(int bitmask) noexcept;

        //- True if the character is a punctuation separator (eg, in ISstream).
        //  Since it could also start a number, SUBTRACT is not included as
        //  a separator.
        //
        //  \param c the character to test, passed as int for consistency with
        //      isdigit, isspace etc.
        inline static bool isseparator(int c) noexcept;


    // Member Functions

    // Status

        //- Return the name for the token type
        static word name(const token::tokenType tokType);

        //- Return the name of the current token type
        inline word name() const;

        //- Return the token type
        inline tokenType type() const noexcept;

        //- Change the token type, for similar types.
        //  This can be used to change between string-like variants
        //  (eg, STRING, VARIABLE, etc)
        //  To change types entirely (eg, STRING to DOUBLE),
        //  use the corresponding assignment operator.
        //
        //  \return true if the change was successful or no change was required
        inline bool setType(const tokenType tokType) noexcept;

        //- The line number for the token
        inline label lineNumber() const noexcept;

        //- Change token line number, return old value
        inline label lineNumber(const label lineNum) noexcept;

        //- True if token is not UNDEFINED or ERROR
        inline bool good() const noexcept;

        //- Token is UNDEFINED
        inline bool undefined() const noexcept;

        //- Token is ERROR
        inline bool error() const noexcept;

        //- Token is BOOL
        inline bool isBool() const noexcept;

        //- Token is FLAG
        inline bool isFlag() const noexcept;

        //- Token is PUNCTUATION
        inline bool isPunctuation() const noexcept;

        //- True if token is PUNCTUATION and equal to parameter
        inline bool isPunctuation(const punctuationToken p) const noexcept;

        //- Token is PUNCTUATION and isseparator
        inline bool isSeparator() const noexcept;

        //- Token is LABEL
        inline bool isLabel() const noexcept;

        //- True if token is LABEL and equal to parameter
        inline bool isLabel(const label value) const noexcept;

        //- Token is FLOAT
        inline bool isFloat() const noexcept;

        //- Token is DOUBLE
        inline bool isDouble() const noexcept;

        //- Token is FLOAT or DOUBLE
        inline bool isScalar() const noexcept;

        //- Token is LABEL, FLOAT or DOUBLE
        inline bool isNumber() const noexcept;

        //- Token is word-variant (WORD, DIRECTIVE)
        inline bool isWord() const noexcept;

        //- Token is word-variant and equal to parameter
        inline bool isWord(const std::string& s) const;

        //- Token is DIRECTIVE (word variant)
        inline bool isDirective() const noexcept;

        //- Token is (quoted) STRING (string variant)
        inline bool isQuotedString() const noexcept;

        //- Token is string-variant
        //- (STRING, EXPRESSION, VARIABLE, VERBATIM, CHAR_DATA)
        inline bool isString() const noexcept;

        //- Token is EXPRESSION (string variant)
        inline bool isExpression() const noexcept;

        //- Token is VARIABLE (string variant)
        inline bool isVariable() const noexcept;

        //- Token is VERBATIM string (string variant)
        inline bool isVerbatim() const noexcept;

        //- Token is CHAR_DATA (string variant)
        inline bool isCharData() const noexcept;

        //- Token is word-variant or string-variant
        //- (WORD, DIRECTIVE, STRING, EXPRESSION, VARIABLE, VERBATIM, CHAR_DATA)
        inline bool isStringType() const noexcept;

        //- Token is COMPOUND
        inline bool isCompound() const noexcept;

        //- True if token is COMPOUND and its type() is equal to parameter
        inline bool isCompound(const word& compoundType) const;

        //- If token is COMPOUND and can be dynamic_cast to the \c Type
        //- return a const pointer to the content.
        //  Return nullptr on any failure.
        template<class Type>
        inline const Type* isCompound() const;

        //- True if token is not UNDEFINED or ERROR. Same as good().
        explicit operator bool() const noexcept { return good(); }


    // Access

        //- Return boolean token value.
        //  Report FatalIOError and return false if token is not BOOL or LABEL
        inline bool boolToken() const;

        //- Return flag bitmask value.
        //  Report FatalIOError and return NO_FLAG if token is not FLAG
        inline int flagToken() const;

        //- Return punctuation character.
        //  Report FatalIOError and return \b \\0 if token is not PUNCTUATION
        inline punctuationToken pToken() const;

        //- Return label value.
        //  Report FatalIOError and return \b 0 if token is not LABEL
        inline label labelToken() const;

        //- Return float value.
        //  Report FatalIOError and return \b 0 if token is not FLOAT
        inline float floatToken() const;

        //- Return double value.
        //  Report FatalIOError and return \b 0 if token is not DOUBLE
        inline double doubleToken() const;

        //- Return float or double value.
        //  Report FatalIOError and return \b 0 if token is not a
        //  FLOAT or DOUBLE
        inline scalar scalarToken() const;

        //- Return label, float or double value.
        //  Report FatalIOError and return \b 0 if token is not a
        //  LABEL, FLOAT or DOUBLE
        inline scalar number() const;

        //- Return const reference to the word contents.
        //  Report FatalIOError and return \b "" if token is not a
        //  WORD or DIRECTIVE
        inline const word& wordToken() const;

        //- Return const reference to the string contents.
        //  Report FatalIOError and return \b "" if token is not a
        //  STRING, EXPRESSION, VARIABLE, VERBATIM, CHAR_DATA
        //  or an upcast WORD or DIRECTIVE
        inline const string& stringToken() const;

        //- Const reference to compound token. Fatal if the wrong type.
        inline const compound& compoundToken() const;

        //- Return reference to compound token. Fatal if the wrong type.
        //- No checks for \em released or \em pending states
        inline compound& refCompoundToken();

        //- Return reference to the underlying encapsulated container
        //- (eg, \c List<label> etc) of a compound token.
        //- No checks for \em released or \em pending states
        template<class Type>
        inline Type& refCompoundToken();

        //- Default construct the specified compound type and read from stream.
        //  A no-op and returns false if compound type is unknown.
        //  Modify the token and return true on success.
        bool readCompoundToken
        (
            const word& compoundType,
            Istream& is,
            const bool readContent = true
        );

        //- Return reference to compound and mark internally as \em released.
        //  Optional Istream parameter only as reference for errors messages.
        compound& transferCompoundToken(const Istream* is = nullptr);

        //- Return reference to compound and mark internally as \em released.
        //  The Istream is used for reference error messages only.
        inline compound& transferCompoundToken(const Istream& is);

        //- Mark the compound as \em released and return a reference
        //- to the underlying encapsulated container - eg, \c List<label>
        //  The Istream is used for reference error messages only.
        template<class Type>
        inline Type& transferCompoundToken(const Istream& is);


    // Edit

        //- Reset token to UNDEFINED and clear any allocated storage
        inline void reset();

        //- Clear token and set to be ERROR.
        inline void setBad();

        //- Swap token contents: type, data, line-number
        inline void swap(token& tok) noexcept;


    // Info

        //- Return info proxy,
        //- for printing token information to a stream
        InfoProxy<token> info() const noexcept { return *this; }


    // Member Operators

    // Assignment

        //- Copy assign
        inline void operator=(const token& tok);

        //- Move assign
        inline void operator=(token&& tok);

        //- Copy assign from punctuation
        inline void operator=(const punctuationToken p);

        //- Copy assign from label
        inline void operator=(const label val);

        //- Copy assign from float
        inline void operator=(const float val);

        //- Copy assign from double
        inline void operator=(const double val);

        //- Copy assign from word content
        inline void operator=(const word& w);

        //- Copy assign from string content
        inline void operator=(const string& str);

        //- Move assign from word content
        inline void operator=(word&& w);

        //- Move assign from string content
        inline void operator=(string&& str);

        //- Assign compound with reference counting to token
        inline void operator=(token::compound* ptr);

        //- Move assign from compound pointer
        inline void operator=(autoPtr<token::compound>&& ptr);


    // Equality

        inline bool operator==(const token& tok) const;
        inline bool operator==(const punctuationToken p) const noexcept;
        inline bool operator==(const label val) const noexcept;
        inline bool operator==(const float val) const noexcept;
        inline bool operator==(const double val) const noexcept;
        inline bool operator==(const std::string& s) const;


    // Inequality

        inline bool operator!=(const token& tok) const;
        inline bool operator!=(const punctuationToken p) const noexcept;
        inline bool operator!=(const label val) const noexcept;
        inline bool operator!=(const float val) const noexcept;
        inline bool operator!=(const double val) const noexcept;
        inline bool operator!=(const std::string& s) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const token& tok);

        friend Ostream& operator<<(Ostream& os, const punctuationToken& pt);
        friend ostream& operator<<(ostream& os, const punctuationToken& pt);

        friend ostream& operator<<(ostream& os, const InfoProxy<token>& ip);


    // Housekeeping

        //- Write access for the token line number
        //  \deprecated(2021-03) - use lineNumber(label)
        label& lineNumber() noexcept { return line_; }

        //- Token is FLOAT
        //  \deprecated(2020-01) - isFloat()
        bool isFloatScalar() const { return isFloat(); };

        //- Token is DOUBLE
        //  \deprecated(2020-01) - isDouble()
        bool isDoubleScalar() const { return isDouble(); }

        //- Return float value.
        //  \deprecated(2020-01) - floatToken()
        float floatScalarToken() const { return floatToken(); }

        //- Return double value.
        //  \deprecated(2020-01) - doubleToken()
        double doubleScalarToken() const { return doubleToken(); }

        //- Deprecated(2017-11) transfer word pointer to the token
        //  \deprecated(2017-11) - use move assign from word
        void operator=(word*) = delete;

        //- Deprecated(2017-11) transfer string pointer to the token
        //  \deprecated(2017-11) - use move assign from string
        void operator=(string*) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IOstream Operators

Istream& operator>>(Istream& is, token& tok);
Ostream& operator<<(Ostream& os, const token::punctuationToken& pt);
ostream& operator<<(ostream& os, const token::punctuationToken& pt);
Ostream& operator<<(Ostream& os, const token::compound& ct);

ostream& operator<<(ostream& os, const InfoProxy<token>& ip);

template<>
Ostream& operator<<(Ostream& os, const InfoProxy<token>& ip);


// Handling of compound types

//- Define compound using \a Type for its name
#define defineCompoundTypeName(Type, UnusedTag)                                \
    defineTemplateTypeNameWithName(token::Compound<Type>, #Type);

//- Define compound using \a Name for its name
#define defineNamedCompoundTypeName(Type, Name)                                \
    defineTemplateTypeNameWithName(token::Compound<Type>, #Name);

//- Add compound to selection tables, lookup using typeName
#define addCompoundToRunTimeSelectionTable(Type, Tag)                          \
    token::compound::addemptyConstructorToTable<token::Compound<Type>>         \
        add##Tag##emptyConstructorToTable_;

//- Add compound to selection tables, lookup as \a Name
#define addNamedCompoundToRunTimeSelectionTable(Type, Tag, Name)               \
    token::compound::addemptyConstructorToTable<token::Compound<Type>>         \
        add##Tag##emptyConstructorToTable_(#Name);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tokenI.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
