/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    vtkPVFoamReader

Description
    reads a dataset in OpenFOAM format

    vtkPVblockMeshReader creates an multiblock dataset.
    It uses the OpenFOAM infrastructure (fvMesh, etc) to handle mesh and
    field data.

SourceFiles
    vtkPVblockMeshReader.cxx

\*---------------------------------------------------------------------------*/
#ifndef vtkPVFoamReader_h
#define vtkPVFoamReader_h

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;
template<class T> class vtkSmartPointer;

// OpenFOAM forward declarations
namespace Foam
{
    class vtkPVFoam;
}


/*---------------------------------------------------------------------------*\
                       Class vtkPVFoamReader Declaration
\*---------------------------------------------------------------------------*/

class vtkPVFoamReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeMacro(vtkPVFoamReader, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream&, vtkIndent) VTK_OVERRIDE;

    static vtkPVFoamReader* New();

    // Description:
    // Get the current timestep and the timestep range.
    vtkGetVector2Macro(TimeStepRange, int);

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // Print some case information
    virtual void PrintInfo();

    // Description:
    // OpenFOAM mesh caching control
    vtkSetMacro(CacheMesh, bool);
    vtkGetMacro(CacheMesh, bool);

    // Description:
    // OpenFOAM refresh times/fields
    virtual void Refresh();

    // Description:
    // OpenFOAM skip/include the 0/ time directory
    vtkSetMacro(SkipZeroTime, bool);
    vtkGetMacro(SkipZeroTime, bool);

    // Description:
    // GUI update control
    vtkSetMacro(UpdateGUI, bool);
    vtkGetMacro(UpdateGUI, bool);

    // Description:
    // OpenFOAM extrapolate internal values onto the patches
    vtkSetMacro(ExtrapolatePatches, bool);
    vtkGetMacro(ExtrapolatePatches, bool);

    // Description:
    // OpenFOAM use vtkPolyhedron instead of decomposing polyhedra
    vtkSetMacro(UseVTKPolyhedron, bool);
    vtkGetMacro(UseVTKPolyhedron, bool);

    // Description:
    // OpenFOAM read sets control
    virtual void SetIncludeSets(bool);
    vtkGetMacro(IncludeSets, bool);

    // Description:
    // OpenFOAM read zones control
    virtual void SetIncludeZones(bool);
    vtkGetMacro(IncludeZones, bool);

    // Description:
    // OpenFOAM display patch names control
    virtual void SetShowPatchNames(bool);
    vtkGetMacro(ShowPatchNames, bool);

    // Description:
    // OpenFOAM display patchGroups
    virtual void SetShowGroupsOnly(bool);
    vtkGetMacro(ShowGroupsOnly, bool);

    // Description:
    // OpenFOAM volField interpolation
    vtkSetMacro(InterpolateVolFields, bool);
    vtkGetMacro(InterpolateVolFields, bool);

    // Description:
    // Get the current timestep
    int  GetTimeStep();

    // Description:
    // Parts selection list control
    virtual vtkDataArraySelection* GetPartSelection();
    int  GetNumberOfPartArrays();
    int  GetPartArrayStatus(const char* name);
    void SetPartArrayStatus(const char* name, int status);
    const char* GetPartArrayName(int index);

    // Description:
    // volField selection list control
    virtual vtkDataArraySelection* GetVolFieldSelection();
    int  GetNumberOfVolFieldArrays();
    int  GetVolFieldArrayStatus(const char* name);
    void SetVolFieldArrayStatus(const char* name, int status);
    const char* GetVolFieldArrayName(int index);

    // Description:
    // pointField selection list control
    virtual vtkDataArraySelection* GetPointFieldSelection();
    int  GetNumberOfPointFieldArrays();
    int  GetPointFieldArrayStatus(const char* name);
    void SetPointFieldArrayStatus(const char* name, int status);
    const char* GetPointFieldArrayName(int index);

    // Description:
    // lagrangianField selection list control
    virtual vtkDataArraySelection* GetLagrangianFieldSelection();
    int  GetNumberOfLagrangianFieldArrays();
    int  GetLagrangianFieldArrayStatus(const char* name);
    void SetLagrangianFieldArrayStatus(const char* name, int status);
    const char* GetLagrangianFieldArrayName(int index);

    // Description:
    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );

    void SelectionModified();


protected:

    //- Construct null
    vtkPVFoamReader();

    //- Destructor
    ~vtkPVFoamReader();

    //- Return information about mesh, times, etc without loading anything
    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    ) VTK_OVERRIDE;

    //- Get the mesh/fields for a particular time
    virtual int RequestData
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    ) VTK_OVERRIDE;

    //- Fill in additional port information
    virtual int FillOutputPortInformation(int, vtkInformation*) VTK_OVERRIDE;

    //- The observer to modify this object when array selections are modified
    vtkCallbackCommand* SelectionObserver;

    //- The file name for this case
    char* FileName;


private:

    //- Disallow default bitwise copy construct
    vtkPVFoamReader(const vtkPVFoamReader&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const vtkPVFoamReader&) = delete;

    //- Add/remove patch names to/from the view
    void updatePatchNamesView(const bool show);

    int TimeStepRange[2];
    bool CacheMesh;
    bool SkipZeroTime;

    bool ExtrapolatePatches;
    bool UseVTKPolyhedron;
    bool IncludeSets;
    bool IncludeZones;
    bool ShowPatchNames;
    bool ShowGroupsOnly;
    bool InterpolateVolFields;

    //- Dummy variable/switch to invoke a reader update
    bool UpdateGUI;

    vtkDataArraySelection* PartSelection;
    vtkDataArraySelection* VolFieldSelection;
    vtkDataArraySelection* PointFieldSelection;
    vtkDataArraySelection* LagrangianFieldSelection;

    //- Cached data for output port0 (experimental!)
    vtkMultiBlockDataSet* output0_;

    //- Backend portion of the reader
    Foam::vtkPVFoam* backend_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
