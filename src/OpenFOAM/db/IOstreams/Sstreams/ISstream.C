/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "ISstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Truncate error message for readability
static const unsigned errLen = 80;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Convert a single character to a word with length 1
inline static Foam::word charToWord(char c)
{
    return Foam::word(std::string(1, c), false);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

char Foam::ISstream::nextValid()
{
    char c = 0;

    while (true)
    {
        // Get next non-whitespace character
        while (get(c) && isspace(c))
        {}

        // Return if stream is bad - ie, previous get() failed
        if (bad() || isspace(c))
        {
            return 0;
        }

        // Is this the start of a C/C++ comment?
        if (c == '/')
        {
            if (!get(c))
            {
                // Cannot get another character - return this one
                return '/';
            }

            if (c == '/')
            {
                // C++ style single-line comment - skip through past end-of-line
                while (get(c) && c != '\n')
                {}
            }
            else if (c == '*')
            {
                // Within a C-style comment
                while (true)
                {
                    // Search for end of C-style comment - '*/'
                    if (get(c) && c == '*')
                    {
                        if (get(c))
                        {
                            if (c == '/')
                            {
                                // matched '*/'
                                break;
                            }
                            else if (c == '*')
                            {
                                // check again
                                putback(c);
                            }
                        }
                    }

                    if (!good())
                    {
                        return 0;
                    }
                }
            }
            else
            {
                // The '/' did not start a C/C++ comment - return it
                putback(c);
                return '/';
            }
        }
        else
        {
            // A valid character - return it
            return c;
        }
    }

    return 0;
}


void Foam::ISstream::readWordToken(token& t)
{
    word val;
    if (read(val).bad())
    {
        t.setBad();
    }
    else if (token::compound::isCompound(val))
    {
        t = token::compound::New(val, *this).ptr();
    }
    else
    {
        t = std::move(val); // Move contents to token
    }
}


Foam::Istream& Foam::ISstream::read(token& t)
{
    static const unsigned maxLen = 128;  // When parsing labels or scalars
    static char buf[maxLen];

    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    // Assume that the streams supplied are in working order.
    // Lines are counted by '\n'

    // Get next 'valid character': i.e. proceed through any whitespace
    // and/or comments until a semantically valid character is found

    char c = nextValid();

    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // Return on error
    if (!c)
    {
        t.setBad();
        return *this;
    }

    // Analyse input starting with this character.
    switch (c)
    {
        // Check for punctuation first - same as token::isSeparator

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
        case token::ADD :
        // NB: token::SUBTRACT handled later as the possible start of a Number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // String: enclosed by double quotes.
        case token::BEGIN_STRING :
        {
            putback(c);

            string val;
            if (read(val).bad())
            {
                t.setBad();
            }
            else
            {
                t = std::move(val); // Move contents to token
            }

            return *this;
        }

        // Possible verbatim string or dictionary functionEntry
        case token::HASH :
        {
            char nextC;
            if (read(nextC).bad())
            {
                // Return lone '#' as word
                t = charToWord(c);
            }
            else if (nextC == token::BEGIN_BLOCK)
            {
                // Verbatim string: #{ ... #}

                string val;
                if (readVerbatim(val).bad())
                {
                    t.setBad();
                }
                else
                {
                    t = std::move(val); // Move contents to token
                    t.setType(token::tokenType::VERBATIMSTRING);
                }
            }
            else
            {
                // Word beginning with '#'. Eg, "#include"
                putback(nextC);
                putback(c);

                readWordToken(t);
            }

            return *this;
        }

        // Dictionary variable (as rvalue)
        case '$':
        {
            char nextC;
            if (read(nextC).bad())
            {
                // Return lone '$' as word
                t = charToWord(c);
            }
            else if (nextC == token::BEGIN_BLOCK)
            {
                // Put back so that "${" is included in the variable
                putback(nextC);
                putback(c);

                string val;
                if (readVariable(val).bad())
                {
                    t.setBad();
                }
                else
                {
                    t = std::move(val); // Move contents to token
                    t.setType(token::tokenType::VARIABLE);
                }
            }
            else
            {
                // Word/variable beginning with '$', but without "{}"

                putback(nextC);
                putback(c);
                readWordToken(t);
            }

            return *this;
        }

        // Number: integer or floating point
        //
        // ideally match the equivalent of this regular expression
        //
        //    /[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([Ee][-+]?[0-9]+)?/
        //
        case '-' :
        case '.' :
        case '0' : case '1' : case '2' : case '3' : case '4' :
        case '5' : case '6' : case '7' : case '8' : case '9' :
        {
            label labelVal = (c != '.'); // used as bool here

            unsigned nChar = 0;
            buf[nChar++] = c;

            // get everything that could resemble a number and let
            // readScalar determine the validity
            while
            (
                is_.get(c)
             && (
                    isdigit(c)
                 || c == '+'
                 || c == '-'
                 || c == '.'
                 || c == 'E'
                 || c == 'e'
                )
            )
            {
                if (labelVal)
                {
                    labelVal = isdigit(c);
                }

                buf[nChar++] = c;
                if (nChar == maxLen)
                {
                    // runaway argument - avoid buffer overflow
                    buf[maxLen-1] = '\0';

                    FatalIOErrorInFunction(*this)
                        << "number '" << buf << "...'\n"
                        << "    is too long (max. " << maxLen << " characters)"
                        << exit(FatalIOError);

                    t.setBad();
                    return *this;
                }
            }
            buf[nChar] = '\0';

            setState(is_.rdstate());
            if (is_.bad())
            {
                t.setBad();
            }
            else
            {
                is_.putback(c);

                if (nChar == 1 && buf[0] == '-')
                {
                    // A single '-' is punctuation
                    t = token::punctuationToken(token::SUBTRACT);
                }
                else if (labelVal && Foam::read(buf, labelVal))
                {
                    t = labelVal;
                }
                else
                {
                    scalar scalarVal;

                    if (readScalar(buf, scalarVal))
                    {
                        // A scalar or too big to fit as a label
                        t = scalarVal;
                    }
                    else
                    {
                        t.setBad();
                    }
                }
            }

            return *this;
        }


        // Should be a word (which can also be a single character)
        default:
        {
            putback(c);
            readWordToken(t);

            return *this;
        }
    }
}


Foam::Istream& Foam::ISstream::read(char& c)
{
    c = nextValid();
    return *this;
}


Foam::Istream& Foam::ISstream::read(word& str)
{
    static const unsigned maxLen = 1024;
    static char buf[maxLen];

    unsigned nChar = 0;
    unsigned depth = 0;  // Track depth of "()" nesting
    char c;

    while (get(c) && word::valid(c))
    {
        if (c == token::BEGIN_LIST)
        {
            ++depth;
        }
        else if (c == token::END_LIST)
        {
            if (depth)
            {
                --depth;
            }
            else
            {
                // Had ')' without a previous '(' ... stop
                break;
            }
        }

        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            buf[errLen] = '\0';

            FatalIOErrorInFunction(*this)
                << "word '" << buf << "...'\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }
    }

    // Terminate string with nul char
    buf[nChar] = '\0';

    // We could probably skip this check
    if (bad())
    {
        buf[errLen] = '\0';

        FatalIOErrorInFunction(*this)
            << "problem while reading word '" << buf << "...' after "
            << nChar << " characters\n"
            << exit(FatalIOError);

        return *this;
    }

    if (nChar == 0)
    {
        FatalIOErrorInFunction(*this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }
    else if (depth)
    {
        IOWarningInFunction(*this)
            << "Missing " << depth << " closing ')' while parsing" << nl << nl
            << buf << nl << endl;
    }

    // Finalize: content already validated, assign without additional checks.
    str.assign(buf, nChar);
    putback(c);

    return *this;
}


Foam::Istream& Foam::ISstream::read(string& str)
{
    static const unsigned maxLen = 1024;
    static char buf[maxLen];

    char c;

    if (!get(c))
    {
        FatalIOErrorInFunction(*this)
            << "cannot read start of string"
            << exit(FatalIOError);

        return *this;
    }

    // Note, we could also handle single-quoted strings here (if desired)
    if (c != token::BEGIN_STRING)
    {
        FatalIOErrorInFunction(*this)
            << "Incorrect start of string character found : " << c
            << exit(FatalIOError);

        return *this;
    }

    unsigned nChar = 0;
    bool escaped = false;

    while (get(c))
    {
        if (c == token::END_STRING)
        {
            if (escaped)
            {
                escaped = false;
                --nChar;    // Overwrite backslash
            }
            else
            {
                // Done reading
                buf[nChar] = '\0';
                str = buf;
                return *this;
            }
        }
        else if (c == token::NL)
        {
            if (escaped)
            {
                escaped = false;
                --nChar;    // Overwrite backslash
            }
            else
            {
                buf[errLen] = buf[nChar] = '\0';

                FatalIOErrorInFunction(*this)
                    << "found '\\n' while reading string \""
                    << buf << "...\""
                    << exit(FatalIOError);

                return *this;
            }
        }
        else if (c == '\\')
        {
            escaped = !escaped;    // toggle state (retains backslashes)
        }
        else
        {
            escaped = false;
        }

        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            buf[errLen] = '\0';

            FatalIOErrorInFunction(*this)
                << "string \"" << buf << "...\"\n"
                << "    is too long (max. " << maxLen << " characters)"
                << exit(FatalIOError);

            return *this;
        }
    }


    // Don't worry about a dangling backslash if string terminated prematurely
    buf[errLen] = buf[nChar] = '\0';

    FatalIOErrorInFunction(*this)
        << "problem while reading string \"" << buf << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::Istream& Foam::ISstream::readVariable(string& str)
{
    static const unsigned maxLen = 1024;
    static char buf[maxLen];

    unsigned nChar = 0;
    unsigned depth = 0;  // Track depth of "{}" nesting
    char c;

    if (!get(c) || c != '$')
    {
        FatalIOErrorInFunction(*this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    buf[nChar++] = c;

    // Read next character to see if '{'
    if (get(c) && c == token::BEGIN_BLOCK)
    {
        buf[nChar++] = c;
        ++depth;  // Starts with '{'

        // Also allow '/' between ${...} blocks for slash-scoping of entries
        while
        (
            get(c)
         && (
                c == token::BEGIN_BLOCK
             || c == token::END_BLOCK
             || word::valid(c) || c == '/'
            )
        )
        {
            if (c == token::BEGIN_BLOCK)
            {
                ++depth;
            }
            else if (c == token::END_BLOCK)
            {
                if (depth)
                {
                    --depth;
                }
                else
                {
                    // Had '}' without a previous '{' ... stop
                    break;
                }
            }

            buf[nChar++] = c;
            if (nChar == maxLen)
            {
                buf[errLen] = '\0';

                FatalIOErrorInFunction(*this)
                    << "variable '" << buf << "...'\n"
                    << "    is too long (max. " << maxLen << " characters)"
                    << exit(FatalIOError);

                return *this;
            }
        }
    }
    else
    {
        buf[nChar++] = c;

        while (get(c) && word::valid(c))
        {
            buf[nChar++] = c;
            if (nChar == maxLen)
            {
                buf[errLen] = '\0';

                FatalIOErrorInFunction(*this)
                    << "variable '" << buf << "...'\n"
                    << "    is too long (max. " << maxLen << " characters)"
                    << exit(FatalIOError);

                return *this;
            }
        }
    }

    // Terminate string with nul char
    buf[nChar] = '\0';

    // we could probably skip this check
    if (bad())
    {
        buf[errLen] = '\0';

        FatalIOErrorInFunction(*this)
            << "problem while reading string '" << buf << "...' after "
            << nChar << " characters\n"
            << exit(FatalIOError);

        return *this;
    }

    if (nChar == 0)
    {
        FatalIOErrorInFunction(*this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }
    else if (depth)
    {
        IOWarningInFunction(*this)
            << "Missing " << depth << " closing '}' while parsing" << nl << nl
            << buf << nl << endl;
    }

    // Finalize
    str = buf;
    putback(c);

    return *this;
}


Foam::Istream& Foam::ISstream::readVerbatim(string& str)
{
    static const unsigned maxLen = 8000;
    static char buf[maxLen];

    unsigned nChar = 0;
    char c;

    str.clear();
    while (get(c))
    {
        if (c == token::HASH)
        {
            char nextC;
            get(nextC);
            if (nextC == token::END_BLOCK)
            {
                // The closing "#}" found
                str.append(buf, nChar);
                return *this;
            }
            else
            {
                putback(nextC);
            }
        }

        buf[nChar++] = c;
        if (nChar == maxLen)
        {
            str.append(buf, nChar);
            nChar = 0;
        }
    }


    // Don't worry about a dangling backslash if string terminated prematurely
    buf[errLen] = buf[nChar] = '\0';

    FatalIOErrorInFunction(*this)
        << "problem while reading string \"" << buf << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::Istream& Foam::ISstream::read(label& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(floatScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


Foam::Istream& Foam::ISstream::read(doubleScalar& val)
{
    is_ >> val;
    setState(is_.rdstate());
    return *this;
}


// read binary block
Foam::Istream& Foam::ISstream::read(char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorInFunction(*this)
            << "stream format not binary"
            << exit(FatalIOError);
    }

    readBegin("binaryBlock");
    is_.read(buf, count);
    readEnd("binaryBlock");

    setState(is_.rdstate());

    return *this;
}


void Foam::ISstream::rewind()
{
    lineNumber_ = 1;      // Reset line number

    stdStream().clear();  // Clear the iostate error state flags
    setGood();            // Sync local copy of iostate

    // pubseekpos() rather than seekg() so that it works with gzstream
    stdStream().rdbuf()->pubseekpos(0, std::ios_base::in);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


std::ios_base::fmtflags Foam::ISstream::flags() const
{
    return is_.flags();
}


std::ios_base::fmtflags Foam::ISstream::flags(const ios_base::fmtflags f)
{
    return is_.flags(f);
}


// ************************************************************************* //
