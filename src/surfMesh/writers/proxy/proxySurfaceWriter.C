/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2009-2011 OpenCFD Ltd.
    Copyright (C) 2011 OpenFOAM Foundation
    Modified code Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "proxySurfaceWriter.H"
#include "MeshedSurfaceProxy.H"
#include "OSspecific.H"
#include "surfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(proxyWriter);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::proxyWriter::proxyWriter(const word& fileExt)
:
    surfaceWriter(),
    fileExtension_(fileExt)
{}


Foam::surfaceWriters::proxyWriter::proxyWriter
(
    const word& fileExt,
    const dictionary& options
)
:
    surfaceWriter(options),
    fileExtension_(fileExt),
    options_(options)
{}


Foam::surfaceWriters::proxyWriter::proxyWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    proxyWriter(outputPath.ext(), options)
{
    surfaceWriter::open(surf, outputPath, parallel);
}


Foam::surfaceWriters::proxyWriter::proxyWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    proxyWriter(outputPath.ext(), options)
{
    surfaceWriter::open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::proxyWriter::write()
{
    checkOpen();

    // Avoid bad values
    if (fileExtension_.empty())
    {
        return fileName::null;
    }

    // Geometry:  rootdir/<TIME>/surfaceName.{extension}

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext(fileExtension_);

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        MeshedSurfaceProxy<face>(surf.points(), surf.faces()).write
        (
            outputFile,
            fileExtension_,
            options_
        );
    }

    wroteGeom_ = true;
    return outputFile;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::proxyWriter);


// ************************************************************************* //
