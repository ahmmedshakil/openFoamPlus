/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::surfMesh

Description
    A surface mesh consisting of general polygon faces.

SourceFiles
    surfMesh.C
    surfMeshClear.C
    surfMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef surfMesh_H
#define surfMesh_H

#include "surfaceRegistry.H"
#include "MeshedSurfaceIOAllocator.H"
#include "PrimitivePatch.H"
#include "SubField.H"
#include "surfZoneIOList.H"

#include "surfFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dimensionSet;
class meshedSurf;
template<class Face> class MeshedSurface;

/*---------------------------------------------------------------------------*\
                          Class surfMesh Declaration
\*---------------------------------------------------------------------------*/

class surfMesh
:
    public  surfaceRegistry,
    private Detail::MeshedSurfaceIOAllocator,
    public  PrimitivePatch<face, ::Foam::UList, ::Foam::SubField<point>, point>
{
public:

    // Public Data Types

        //- Enumeration defining the state of the mesh after a read update.
        //  Used for post-processing applications, where the mesh
        //  needs to update based on the files written in time directories
        enum readUpdateState
        {
            UNCHANGED,
            POINTS_MOVED,
            TOPO_CHANGE,
            TOPO_PATCH_CHANGE
        };


private:

    // Private Typedefs

        typedef Detail::MeshedSurfaceIOAllocator Allocator;

        typedef PrimitivePatch
        <
            face,
            ::Foam::UList,
            ::Foam::SubField<point>,
            point
        >
        MeshReference;


    // Demand driven private data

        //- The surface zones
        surfZoneIOList surfZones_;


    // Private Member Functions

        //- No copy construct
        surfMesh(const surfMesh&) = delete;

        //- No copy assignment
        void operator=(const surfMesh&) = delete;


protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return Allocator::storedIOPoints();
        }

        //- Non-const access to the faces
        faceList& storedFaces()
        {
            return Allocator::storedIOFaces();
        }

        //- Non-const access to the surface zones
        surfZoneList& storedZones()
        {
            return surfZones_;
        }

        //- Update references to storedFaces
        virtual void updateFacesRef();

        //- Update references to storedPoints
        virtual void updatePointsRef();

        //- Update references to storedPoints/storedFaces
        virtual void updateRefs();


public:

    // Public typedefs

    //- Typedef required for GeoMesh
    typedef surfMesh Mesh;

    //- Placeholder only, but do not remove - it is needed for GeoMesh
    typedef bool BoundaryMesh;


    //- Runtime type information
    TypeName("surfMesh");

    //- Return the mesh sub-directory name (normally "surfMesh")
    static word meshSubDir;


    // Constructors

        //- Read construct from IOobject.
        //  Writing = NO_WRITE
        explicit surfMesh(const IOobject& io);

        //- Read construct from IOobject, with alternative surface name
        //  Writing = NO_WRITE
        surfMesh(const IOobject& io, const word& surfName);

        //- Construct null with specified name on the given registry.
        surfMesh(const word& surfName, const objectRegistry& obr);

        //- Copy construct from MeshedSurface<face>
        surfMesh
        (
            const IOobject& io,
            const MeshedSurface<face>& surf,
            const word& surfName = word::null
        );

        //- Move construct from MeshedSurface<face>
        surfMesh
        (
            const IOobject& io,
            MeshedSurface<face>&& surf,
            const word& surfName = word::null
        );


    //- Destructor
    virtual ~surfMesh();


    // Member Functions

    // Database

        //- Return the local mesh directory (dbDir()/meshSubDir)
        fileName meshDir() const;

        //- Return the current instance directory for points
        //  Used when constructing geometric mesh data dependent on points
        const fileName& pointsInstance() const;

        //- Return the current instance directory for faces
        const fileName& facesInstance() const;

        //- Set the instance for mesh files
        void setInstance
        (
            const fileName& inst,
            IOobject::writeOption wOpt = IOobject::AUTO_WRITE
        );

        //- Adjust the write option for all components
        void setWriteOption(IOobject::writeOption wOpt);


    // Access

        //- Return the number of raw points
        virtual label nPoints() const;

        //- Return the number of raw faces
        virtual label nFaces() const;

        //- Return number of faces
        virtual label size() const
        {
            return nFaces();
        }


        //- Return points
        virtual const pointField& points() const;

        //- Return faces
        virtual const faceList& faces() const;

        //- Return surface zones
        virtual const surfZoneList& surfZones() const
        {
            return surfZones_;
        }

        //- Return face area vectors (normals)
        inline const vectorField& Sf() const
        {
            return MeshReference::faceAreas();
        }

        //- Return face area magnitudes
        inline const scalarField& magSf() const
        {
            return MeshReference::magFaceAreas();
        }

        //- Face centres
        inline const vectorField& Cf() const
        {
            return MeshReference::faceCentres();
        }


    // Zones

        //- Add surface zones, optionally validating the zone coverage
        void addZones
        (
            const surfZoneList& zones,
            bool validate = true
        );

        //- Remove surface zones
        void removeZones();

        //- Check the surface zone definitions
        void checkZones();


    // Modification

        //- Update with new contents
        void copySurface
        (
            const pointField& points,
            const faceList& faces,
            bool validate=false
        );

        //- Update with new contents
        void copySurface
        (
            const meshedSurf& surf,
            bool validate=false
        );

        //- Update with new contents
        void copySurface
        (
            const MeshedSurface<face>& surf,
            bool validate=false
        );

        //- Transfer the contents of the argument and annul the argument
        //  Optionally validate the zone coverage.
        void transfer(MeshedSurface<face>& surf, bool validate=false);

        //- Update mesh based on the files saved in time directories
        virtual readUpdateState readUpdate();


    // Fields

        //- Copy/store named field as face or point data (template parameter).
        //
        //  Default is face-data (surfGeoMesh as template).
        template<class Type, class GeoMeshType = surfGeoMesh>
        void storeField
        (
            const word& fieldName,
            const dimensionSet& dims,
            const Field<Type>& values
        );

        //- Move/store named field as face or point data (template parameter).
        //
        //  Default is face-data (surfGeoMesh as template).
        template<class Type, class GeoMeshType = surfGeoMesh>
        void storeField
        (
            const word& fieldName,
            const dimensionSet& dims,
            Field<Type>&& values
        );


    // Writing

        //- Avoid masking the normal objectRegistry write
        using surfaceRegistry::write;

        //- Write all components using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool valid
        ) const;


        //- Write to file, choosing writer based on its extension.
        //  Uses MeshedSurfaceProxy for writing.
        void write
        (
            const fileName& name,
            const dictionary& options = dictionary::null
        ) const;

        //- Write to file, choosing writer given extension.
        //  Uses MeshedSurfaceProxy for writing.
        void write
        (
            const fileName& name,
            const word& ext,
            const dictionary& options = dictionary::null
        ) const;


    // Storage management

        //- Release the geometry and return as a MeshedSurface<face>.
        autoPtr<MeshedSurface<face>> releaseGeom();

        //- Clear geometry
        void clearGeom();

        //- Clear addressing
        void clearAddressing();

        //- Clear all geometry and addressing unnecessary for CFD
        void clearOut();

        //- Clear primitive data (points, faces and cells)
        void clearPrimitives();

        //- Clear stored fields
        void clearFields();

        //- Remove all files from mesh instance
        void removeFiles(const fileName& instanceDir) const;

        //- Remove all files from mesh instance()
        void removeFiles() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
