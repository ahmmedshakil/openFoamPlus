/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Application
    blockMesh

Group
    grpMeshGenerationUtilities

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
      - \c system/blockMeshDict
      - \c system/\<region\>/blockMeshDict
      - \c constant/polyMesh/blockMeshDict
      - \c constant/\<region\>/polyMesh/blockMeshDict

Usage
    \b blockMesh [OPTION]

    Options:
      - \par -blockTopology
        Write the topology as a set of edges in OBJ format.

      - \par -region \<name\>
        Specify an alternative mesh region.

      - \par -dict \<filename\>
        Specify alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "blockTopology",
        "write block edges and centres as .obj files"
    );
    argList::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the blockMesh description"
    );

    argList::addNote
    (
        "Block description\n"
        "\n"
        "  For a given block, the correspondence between the ordering of\n"
        "  vertex labels and face labels is shown below.\n"
        "  For vertex numbering in the sequence 0 to 7 (block, centre):\n"
        "    faces 0 (f0) and 1 are left and right, respectively;\n"
        "    faces 2 and 3 are bottom and top;\n"
        "    and faces 4 and 5 are front the back:\n"
        "\n"
        "           4 ---- 5\n"
        "      f3   |\\     |\\   f5\n"
        "      |    | 7 ---- 6   \\\n"
        "      |    0 |--- 1 |    \\\n"
        "      |     \\|     \\|    f4\n"
        "      f2     3 ---- 2\n"
        "\n"
        "            f0 ----- f1\n"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("blockMeshDict");

    word regionName;
    word regionPath;

    // Check if the region is specified otherwise mesh the default region
    if (args.optionReadIfPresent("region", regionName, polyMesh::defaultRegion))
    {
        Info<< nl << "Generating mesh for region " << regionName << endl;
        regionPath = regionName;
    }

    // Search for the appropriate blockMesh dictionary....

    fileName dictPath;

    // Check if the dictionary is specified on the command-line
    if (args.optionReadIfPresent("dict", dictPath))
    {
        if (isDir(dictPath))
        {
            dictPath = dictPath / dictName;
        }
    }
    // Check if dictionary is present in the constant directory
    else if
    (
        exists
        (
            runTime.path()/runTime.constant()
           /regionPath/polyMesh::meshSubDir/dictName
        )
    )
    {
        dictPath =
            runTime.constant()
           /regionPath/polyMesh::meshSubDir/dictName;

        // Warn that constant/polyMesh/blockMesh was selected instead of
        // system/blockMesh
        WarningIn(args[0])
            << "Using the old blockMeshDict location: "
            << dictPath << nl
            << "    instead of the default location:  "
            << runTime.system()/regionPath/dictName << nl
            << endl;
    }
    // Otherwise assume the dictionary is present in the system directory
    else
    {
        dictPath = runTime.system()/regionPath/dictName;
    }

    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << meshDictIO.objectPath()
            << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from\n    "
        << meshDictIO.objectPath() << endl;

    IOdictionary meshDict(meshDictIO);
    blockMesh blocks(meshDict, regionName);


    if (args.optionFound("blockTopology"))
    {
        // Write mesh as edges.
        {
            fileName objMeshFile("blockTopology.obj");

            OFstream str(runTime.path()/objMeshFile);

            Info<< nl << "Dumping block structure as Lightwave obj format"
                << " to " << objMeshFile << endl;

            blocks.writeTopology(str);
        }

        // Write centres of blocks
        {
            fileName objCcFile("blockCentres.obj");

            OFstream str(runTime.path()/objCcFile);

            Info<< nl << "Dumping block centres as Lightwave obj format"
                << " to " << objCcFile << endl;

            const polyMesh& topo = blocks.topology();

            const pointField& cellCentres = topo.cellCentres();

            forAll(cellCentres, celli)
            {
                //point cc = b.blockShape().centre(b.points());
                const point& cc = cellCentres[celli];

                str << "v " << cc.x() << ' ' << cc.y() << ' ' << cc.z() << nl;
            }
        }

        Info<< nl << "end" << endl;

        return 0;
    }


    Info<< nl << "Creating polyMesh from blockMesh" << endl;

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        xferCopy(blocks.points()),           // could we re-use space?
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        defaultFacesName,
        defaultFacesType
    );


    // Read in a list of dictionaries for the merge patch pairs
    if (meshDict.found("mergePatchPairs"))
    {
        List<Pair<word>> mergePatchPairs
        (
            meshDict.lookup("mergePatchPairs")
        );

        #include "mergePatchPairs.H"
    }
    else
    {
        Info<< nl << "There are no merge patch pairs edges" << endl;
    }


    // Set any cellZones (note: cell labelling unaffected by above
    // mergePatchPairs)

    label nZones = blocks.numZonedBlocks();

    if (nZones > 0)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label>> zoneCells(nZones);

        // Running cell counter
        label celli = 0;

        // Largest zone so far
        label freeZoneI = 0;

        forAll(blocks, blockI)
        {
            const block& b = blocks[blockI];
            const List<FixedList<label, 8>> blockCells = b.cells();
            const word& zoneName = b.zoneName();

            if (zoneName.size())
            {
                HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

                label zoneI;

                if (iter == zoneMap.end())
                {
                    zoneI = freeZoneI++;

                    Info<< "    " << zoneI << '\t' << zoneName << endl;

                    zoneMap.insert(zoneName, zoneI);
                }
                else
                {
                    zoneI = iter();
                }

                forAll(blockCells, i)
                {
                    zoneCells[zoneI].append(celli++);
                }
            }
            else
            {
                celli += b.cells().size();
            }
        }


        List<cellZone*> cz(zoneMap.size());

        Info<< nl << "Writing cell zones as cellSets" << endl;

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI] = new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                zoneI,
                mesh.cellZones()
            );

            // Write as cellSet for ease of processing
            cellSet cset(mesh, iter.key(), zoneCells[zoneI].shrink());
            cset.write();
        }

        mesh.pointZones().setSize(0);
        mesh.faceZones().setSize(0);
        mesh.cellZones().setSize(0);
        mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh" << endl;
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


    //
    // write some information
    //
    {
        const polyPatchList& patches = mesh.boundaryMesh();

        Info<< "----------------" << nl
            << "Mesh Information" << nl
            << "----------------" << nl
            << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
            << "  " << "nPoints: " << mesh.nPoints() << nl
            << "  " << "nCells: " << mesh.nCells() << nl
            << "  " << "nFaces: " << mesh.nFaces() << nl
            << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

        Info<< "----------------" << nl
            << "Patches" << nl
            << "----------------" << nl;

        forAll(patches, patchi)
        {
            const polyPatch& p = patches[patchi];

            Info<< "  " << "patch " << patchi
                << " (start: " << p.start()
                << " size: " << p.size()
                << ") name: " << p.name()
                << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //