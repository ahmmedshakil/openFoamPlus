/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "OFFsurfaceFormat.H"
#include "clock.H"
#include "Fstream.H"
#include "StringStream.H"
#include "faceTraits.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OFFsurfaceFormat<Face>::OFFsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OFFsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    // Read header
    string hdr = this->getLineNoComment(is);
    if (hdr != "OFF")
    {
        FatalErrorInFunction
            << "OFF file " << filename << " does not start with 'OFF'"
            << exit(FatalError);
    }


    // get dimensions
    label nPoints, nElems, nEdges;

    string line = this->getLineNoComment(is);
    {
        IStringStream lineStream(line);
        lineStream >> nPoints >> nElems >> nEdges;
    }

    // Read points
    pointField pointLst(nPoints);
    forAll(pointLst, pointi)
    {
        scalar x, y, z;
        line = this->getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream >> x >> y >> z;
        }
        pointLst[pointi] = point(x, y, z);
    }

    // Read faces - ignore optional zone information
    // use a DynamicList for possible on-the-fly triangulation
    DynamicList<Face>  dynFaces(nElems);

    for (label facei = 0; facei < nElems; ++facei)
    {
        line = this->getLineNoComment(is);

        {
            IStringStream lineStream(line);

            label nVerts;
            lineStream >> nVerts;

            List<label> verts(nVerts);

            forAll(verts, vertI)
            {
                lineStream >> verts[vertI];
            }

            const labelUList& f = static_cast<const labelUList&>(verts);

            if (faceTraits<Face>::isTri() && f.size() > 3)
            {
                // simple face triangulation about f[0]
                // cannot use face::triangulation (points may be incomplete)
                for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
                {
                    label fp2 = f.fcIndex(fp1);

                    dynFaces.append(Face{f[0], f[fp1], f[fp2]});
                }
            }
            else
            {
                dynFaces.append(Face(f));
            }
        }
    }

    // Transfer to normal lists, no zone information
    MeshedSurface<Face> surf(std::move(pointLst), std::move(dynFaces));

    this->swap(surf);

    return true;
}


template<class Face>
void Foam::fileFormats::OFFsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf,
    const dictionary&
)
{
    const UList<point>& pointLst = surf.points();
    const UList<Face>&  faceLst  = surf.surfFaces();
    const UList<label>& faceMap  = surf.faceMap();
    const UList<surfZone>& zoneLst = surf.surfZones();

    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing " << filename
            << exit(FatalError);
    }

    // Write header
    os  << "OFF" << endl
        << "# Geomview OFF file written " << clock::dateTime().c_str() << nl
        << nl
        << "# points : " << pointLst.size() << nl
        << "# faces  : " << faceLst.size() << nl
        << "# zones  : " << zoneLst.size() << nl;

    // Print zone names as comment
    forAll(zoneLst, zoneI)
    {
        os  << "#   " << zoneI << "  " << zoneLst[zoneI].name()
            << "  (nFaces: " << zoneLst[zoneI].size() << ")" << nl;
    }

    os  << nl
        << "# nPoints  nFaces  nEdges" << nl
        << pointLst.size() << ' ' << faceLst.size() << ' ' << 0 << nl
        << nl
        << "# <points count=\"" << pointLst.size() << "\">" << endl;

    // Write vertex coords
    forAll(pointLst, ptI)
    {
        os  << pointLst[ptI].x() << ' '
            << pointLst[ptI].y() << ' '
            << pointLst[ptI].z() << " #" << ptI << endl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <faces count=\"" << faceLst.size() << "\">" << endl;

    label faceIndex = 0;
    forAll(zoneLst, zoneI)
    {
        os << "# <zone name=\"" << zoneLst[zoneI].name() << "\">" << endl;

        const label nLocalFaces = zoneLst[zoneI].size();

        if (surf.useFaceMap())
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                os << f.size();
                for (const label verti : f)
                {
                    os << ' ' << verti;
                }

                // add optional zone information
                os << ' ' << zoneI << endl;
            }
        }
        else
        {
            for (label i=0; i<nLocalFaces; ++i)
            {
                const Face& f = faceLst[faceIndex++];

                os << f.size();
                for (const label verti : f)
                {
                    os << ' ' << verti;
                }

                // add optional zone information
                os << ' ' << zoneI << endl;
            }
        }
        os << "# </zone>" << endl;
    }
    os << "# </faces>" << endl;
}


// ************************************************************************* //
