/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "ensightSurfaceReader.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightSurfaceReader, 0);
    addToRunTimeSelectionTable(surfaceReader, ensightSurfaceReader, fileName);
}


void Foam::ensightSurfaceReader::skip(const label n, Istream& is) const
{
    label i = 0;
    token tok;
    while (is.good() && (i < n))
    {
        is >> tok;
        ++i;

        DebugInfo
            << "Skipping token " << tok << nl;
    }

    if (i != n)
    {
        WarningInFunction
            << "Requested to skip " << n << "tokens, but stream exited after "
            << i << " tokens. Last token read: " << tok
            << nl;
    }
}


void Foam::ensightSurfaceReader::readLine(IFstream& is, string& buffer) const
{
    buffer.clear();
    while (is.good() && buffer.empty())
    {
        is.getLine(buffer);
    }
}


void Foam::ensightSurfaceReader::debugSection
(
    const word& expected,
    IFstream& is
) const
{
    string actual;
    readLine(is, actual);

    if (expected != actual)
    {
        FatalIOErrorInFunction(is)
            << "Expected section header '" << expected
            << "' but read the word '" << actual << "'"
            << exit(FatalIOError);
    }

    DebugInfo
        << "Read section header: " << expected << nl;
}


void Foam::ensightSurfaceReader::readGeometryHeader(ensightReadFile& is) const
{
    // Binary flag string if applicable
    is.readBinaryHeader();

    string buffer;

    // Ensight Geometry File
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;

    // Description - 1
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;

    // Node info
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;

    // Element info
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;

    // Part
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;

    // Part number
    label ibuffer;
    is.read(ibuffer);
    DebugInfo<< "ibuffer: " << ibuffer << nl;

    // Description - 2
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;

    // Coordinates
    is.read(buffer);
    DebugInfo<< "buffer: " << buffer << nl;
}


void Foam::ensightSurfaceReader::readCase(IFstream& is)
{
    DebugInFunction<< nl;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    string buffer;

    // Read the file

    debugSection("FORMAT", is);
    readLine(is, buffer);  // type: ensight gold

    debugSection("GEOMETRY", is);
    readLine(is, buffer);

    // GEOMETRY with any of these
    //     model: 1 xxx.0000.mesh
    //     model: xxx.0000.mesh
    //     model: data/directory/geometry
    //
    // - use the last entry
    meshFileName_ = stringOps::splitSpace(buffer).last().str();

    debugSection("VARIABLE", is);

    // Read the field description
    DynamicList<word> fieldNames(10);
    DynamicList<string> fieldFileNames(10);

    while (is.good())
    {
        readLine(is, buffer);

        if (buffer == "TIME")
        {
            break;
        }

        // Read the field name and associated file name. Eg,
        //     scalar per element: 1  p  data/********/p

        const auto parsed = stringOps::splitSpace(buffer);

        if (buffer.find(':') == string::npos || parsed.size() < 4)
        {
            WarningInFunction
                << "Error reading field file name. Current buffer: "
                << buffer << endl;
            continue;
        }
        else if (debug)
        {
            Info<<"variable line: " << parsed.size();
            for (const auto& s : parsed)
            {
                Info<<" " << s.str();
            }
            Info<<nl;
        }

        fieldNames.append(parsed[parsed.size()-2].str());
        fieldFileNames.append(parsed.last().str());
    }
    fieldNames_.transfer(fieldNames);
    fieldFileNames_.transfer(fieldFileNames);

    DebugInfo
        << "fieldNames: " << fieldNames_ << nl
        << "fieldFileNames: " << fieldFileNames_ << nl;

    // Start reading time information
    readLine(is, buffer);                       // time set: <int>

    readLine(is, buffer);
    readFromLine(3, buffer, nTimeSteps_);       // number of steps: <int>
    readLine(is, buffer);
    readFromLine(3, buffer, timeStartIndex_);   // filename start number: <int>
    readLine(is, buffer);
    readFromLine(2, buffer, timeIncrement_);    // filename increment: <int>

    DebugInfo
        << "nTimeSteps: " << nTimeSteps_ << nl
        << "timeStartIndex: " << timeStartIndex_ << nl
        << "timeIncrement: " << timeIncrement_ << nl;

    // Read the time values
    readLine(is, buffer); // time values:
    timeValues_.setSize(nTimeSteps_);
    for (label i = 0; i < nTimeSteps_; i++)
    {
        scalar t(readScalar(is));

        timeValues_[i].value() = t;
        // TODO: use character representation of t directly instead of
        // regenerating from scalar value
        timeValues_[i].name() = Foam::name(t);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightSurfaceReader::ensightSurfaceReader(const fileName& fName)
:
    surfaceReader(fName),
    streamFormat_(IOstream::ASCII),
    baseDir_(fName.path()),
    meshFileName_(),
    fieldNames_(),
    fieldFileNames_(),
    nTimeSteps_(0),
    timeStartIndex_(0),
    timeIncrement_(1),
    timeValues_(),
    surfPtr_(nullptr)
{
    IFstream is(fName);
    readCase(is);
}


// * * * * * * * * * * * * * Public Member Functions   * * * * * * * * * * * //

const Foam::meshedSurface& Foam::ensightSurfaceReader::geometry()
{
    DebugInFunction<< nl;

    if (!surfPtr_.valid())
    {
        IFstream isBinary(baseDir_/meshFileName_, IOstream::BINARY);

        if (!isBinary.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << isBinary.name()
                << exit(FatalError);
        }

        streamFormat_ = IOstream::BINARY;
        {
            istream& is = isBinary.stdStream();

            char buffer[80];
            is.read(buffer, 80);

            char test[80];
            label nChar = 0;
            for (label i = 0; i < 80; ++i)
            {
                if (buffer[i] == '\0')
                {
                    break;
                }
                test[i] = buffer[i];
                nChar++;
            }

            string testStr(test, nChar);

            if
            (
                (testStr.find("binary", 0) == string::npos)
             && (testStr.find("Binary", 0) == string::npos)
            )
            {
                streamFormat_ = IOstream::ASCII;
            }
        }

        if (debug)
        {
            Info<< "stream format: ";
            if (streamFormat_ == IOstream::ASCII)
            {
                Info<< "ascii" << endl;
            }
            else
            {
                Info<< "binary" << endl;
            }
        }


        ensightReadFile is(baseDir_/meshFileName_, streamFormat_);

        DebugInfo
            << "File: " << is.name() << nl;

        readGeometryHeader(is);

        label nPoints;
        is.read(nPoints);

        DebugInfo
            << "nPoints: " << nPoints << nl;

        pointField points(nPoints);
        {
            scalarField x(nPoints);
            for (label dir = 0; dir < 3; dir++)
            {
                forAll(points, pointI)
                {
                    is.read(x[pointI]);
                }

                points.replace(dir, x);
            }
        }


        // Read faces - may be a mix of tris, quads and polys
        DynamicList<face> faces(ceil(nPoints/3));
        DynamicList<Tuple2<string, label>> schema(faces.size());
        string faceType = "";
        label nFace = 0;
        while (is.good()) // (is.peek() != EOF)
        {
            is.read(faceType);

            if (!is.good())
            {
                break;
            }

            DebugInfo
                << "faceType: " << faceType << endl;

            if (faceType == "tria3")
            {
                is.read(nFace);

                const label np = 3;
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    face f(np);
                    for (label fpI = 0; fpI < np; fpI++)
                    {
                        is.read(f[fpI]);
                    }

                    faces.append(f);
                }
            }
            else if (faceType == "quad4")
            {
                is.read(nFace);

                const label np = 4;
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    face f(np);
                    for (label fpI = 0; fpI < np; fpI++)
                    {
                        is.read(f[fpI]);
                    }

                    faces.append(f);
                }
            }
            else if (faceType == "nsided")
            {
                is.read(nFace);

                labelList np(nFace);
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    is.read(np[faceI]);
                }
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    face f(np[faceI]);
                    for (label fpI = 0; fpI < f.size(); ++fpI)
                    {
                        is.read(f[fpI]);
                    }

                    faces.append(f);
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unknown face type: " << faceType
                        << ".  Aborting read and continuing with current "
                        << "elements only" << endl;
                }

                break;
            }
            schema.append(Tuple2<string, label>(faceType, nFace));
        }

        schema_.transfer(schema);

        DebugInfo
            << "read nFaces: " << faces.size() << nl
            << "file schema: " << schema_ << nl;

        // Convert from 1-based Ensight addressing to 0-based OF addressing
        for (face& f : faces)
        {
            for (label& pointi : f)
            {
                --pointi;
            }
        }

        surfPtr_.reset(new meshedSurface(std::move(points), std::move(faces)));
    }

    return *surfPtr_;
}


Foam::instantList Foam::ensightSurfaceReader::times() const
{
    return timeValues_;
}


Foam::wordList Foam::ensightSurfaceReader::fieldNames
(
    const label timeIndex
) const
{
    return fieldNames_;
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const scalar& refValue
) const
{
    return readField<scalar>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const vector& refValue
) const
{
    return readField<vector>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const sphericalTensor& refValue
) const
{
    return readField<sphericalTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const symmTensor& refValue
) const
{
    return readField<symmTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const tensor& refValue
) const
{
    return readField<tensor>(timeIndex, fieldIndex);
}


// ************************************************************************* //
