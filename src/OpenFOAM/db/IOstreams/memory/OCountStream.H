/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::OSCountStream

Description
    An output stream for calculating byte counts.

\*---------------------------------------------------------------------------*/

#ifndef OScountStream_H
#define OScountStream_H

#include "OSstream.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class countstreambuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf class for determining byte counts
class countstreambuf
:
    public std::streambuf
{
    //- The number of bytes counted
    std::streamsize size_;

protected:

    //- Handle output counting via overflow
    virtual int overflow(int c = EOF)
    {
        if (c != EOF)
        {
            ++size_;
        }
        return c;
    }

    //- Set position pointer to absolute position
    //  For the counter, any positioning is ignored and it always acts like
    //  seekpos(0), which resets the count.
    virtual std::streampos seekpos
    (
        std::streampos,
        std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
    )
    {
        size_ = 0;
        return 0;
    }

public:

    //- Construct null, or with precount size
    countstreambuf(std::streamsize precount=0)
    :
        size_(precount)
    {}

    //- \return The number of bytes counted
    std::streamsize size() const
    {
        return size_;
    }

    //- Some information about the number of bytes counted
    inline void printBufInfo(Ostream& os) const
    {
        os << "count=" << size_;
    }
};


/*---------------------------------------------------------------------------*\
                        Class ocountstream Declaration
\*---------------------------------------------------------------------------*/

//- Trivial output stream for calculating byte counts.
//  Since all output values are discarded, it can be used as a /dev/null
//  output buffer as well.
class ocountstream
:
    virtual public std::ios,
    protected countstreambuf,
    public std::ostream
{
public:

    //- Construct null
    ocountstream()
    :
        countstreambuf(),
        std::ostream(static_cast<countstreambuf*>(this))
    {}


    //- \return The number of bytes counted
    using countstreambuf::size;

    //- Rewind the stream, reset the count
    void rewind()
    {
        this->pubseekpos(0, std::ios_base::out);
        clear(); // for safety, clear any old errors
    }
};


namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::OCountStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer allocator for counting
class OCountStreamAllocator
{
protected:

    // Protected data

        typedef std::ostream stream_type;

        //- The stream buffer
        countstreambuf buf_;

        //- The output stream
        stream_type stream_;


    // Constructors

        //- Construct null, or with precount size
        OCountStreamAllocator(std::streamsize precount=0)
        :
            buf_(precount),
            stream_(&buf_)
        {}


    // Protected Member Functions

        void printBufInfo(Ostream& os) const
        {
            buf_.printBufInfo(os);
        }

public:

    // Member Functions

        //- The number of bytes counted
        std::streamsize size() const
        {
            return buf_.size();
        }

        //- Rewind the stream, reset the count
        void rewind()
        {
            buf_.pubseekpos(0);
            stream_.clear(); // for safety, clear any old errors
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                        Class OCountStream Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for calculating byte counts
class OCountStream
:
    public Detail::OCountStreamAllocator,
    public OSstream
{
    typedef Detail::OCountStreamAllocator allocator_type;

public:

    // Constructors

        //- Construct and set stream status
        OCountStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(),
            OSstream(stream_, "output", format, version)
        {}

        //- Copy construct
        OCountStream(const OCountStream& os)
        :
            allocator_type(os.size()),
            OSstream(stream_, os.name(), os.format(), os.version())
        {}


    // Member Functions

        //- Rewind the stream, reset the count, clearing any old errors
        virtual void rewind()
        {
            allocator_type::rewind();
            setGood();  // resynchronize with internal state
        }

        //- Print description to Ostream
        virtual void print(Ostream& os) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
