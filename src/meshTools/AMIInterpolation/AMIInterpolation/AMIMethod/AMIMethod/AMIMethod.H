/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::AMIMethod

Description
    Base class for Arbitrary Mesh Interface (AMI) methods

SourceFiles
    AMIMethod.C

\*---------------------------------------------------------------------------*/

#ifndef AMIMethod_H
#define AMIMethod_H

#include "className.H"
#include "DynamicList.H"
#include "faceAreaIntersect.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"
#include "treeBoundBoxList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class SourcePatch, class TargetPatch> class AMIInterpolation;

/*---------------------------------------------------------------------------*\
                          Class AMIMethod Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class AMIMethod
{

private:

    // Private Member Functions

        //- No copy construct
        AMIMethod(const AMIMethod&) = delete;

        //- No copy assignment
        void operator=(const AMIMethod&) = delete;


protected:

    //- Local typedef to octree tree-type
    typedef treeDataPrimitivePatch<TargetPatch> treeType;


    // Protected data

        //- Reference to source patch
        const SourcePatch& srcPatch_;

        //- Reference to target patch
        const TargetPatch& tgtPatch_;

        //- Flag to indicate that the two patches are co-directional and
        //- that the orientation of the target patch should be reversed
        const bool reverseTarget_;

        //- Flag to indicate that the two patches must be matched/an overlap
        //- exists between them
        const bool requireMatch_;

        //- Source face areas
        List<scalar> srcMagSf_;

        //- Target face areas
        List<scalar> tgtMagSf_;

        //- Labels of faces that are not overlapped by any target faces
        //- (should be empty for correct functioning)
        labelList srcNonOverlap_;

        //- Octree used to find face seeds
        autoPtr<indexedOctree<treeType>> treePtr_;

        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;


    // Protected Member Functions

        // Helper functions

            //- Check AMI patch coupling
            void checkPatches() const;

            //- Initialise and return true if all ok
            bool initialise
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label& srcFacei,
                label& tgtFacei
            );

            //- Write triangle intersection to OBJ file
            void writeIntersectionOBJ
            (
                const scalar area,
                const face& f1,
                const face& f2,
                const pointField& f1Points,
                const pointField& f2Points
            ) const;


        // Common AMI method functions

            //- Reset the octree for the target patch face search
            void resetTree();

            //- Find face on target patch that overlaps source face
            label findTargetFace(const label srcFacei) const;

            //- Add faces neighbouring facei to the ID list
            void appendNbrFaces
            (
                const label facei,
                const TargetPatch& patch,
                const DynamicList<label>& visitedFaces,
                DynamicList<label>& faceIDs
            ) const;

            //- Helper function to decompose a patch
            template<class PatchType>
            void triangulatePatch
            (
                const PatchType& patch,
                List<DynamicList<face>>& tris,
                List<scalar>& magSf
            ) const;


public:

    //- Runtime type information
    TypeName("AMIMethod");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        AMIMethod,
        components,
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget,
            const bool requireMatch
        ),
        (
            srcPatch,
            tgtPatch,
            triMode,
            reverseTarget,
            requireMatch
        )
    );

    //- Construct from components
    AMIMethod
    (
        const SourcePatch& srcPatch,
        const TargetPatch& tgtPatch,
        const faceAreaIntersect::triangulationMode& triMode,
        const bool reverseTarget,
        const bool requireMatch
    );

    //- Selector
    static autoPtr<AMIMethod> New
    (
        const word& methodName,
        const SourcePatch& srcPatch,
        const TargetPatch& tgtPatch,
        const faceAreaIntersect::triangulationMode& triMode,
        const bool reverseTarget,
        const bool requireMatch
    );


    //- Destructor
    virtual ~AMIMethod() = default;


    // Member Functions

        // Access

            //- Labels of faces that are not overlapped by any target faces
            //  Note: this should be empty for correct functioning
            inline const labelList& srcNonOverlap() const;

            //- Flag to indicate that interpolation patches are conformal
            virtual bool conformal() const;

            //- Return const access to source patch face areas
            inline const List<scalar>& srcMagSf() const;

            //- Return access to source patch face areas
            inline List<scalar>& srcMagSf();

            //- Return const access to target patch face areas
            inline const List<scalar>& tgtMagSf() const;

            //- Return access to target patch face areas
            inline List<scalar>& tgtMagSf();


        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            ) = 0;

            //- Set the face areas for parallel runs
            virtual void setMagSf
            (
                const TargetPatch& tgtPatch,
                const mapDistribute& map,
                scalarList& srcMagSf,
                scalarList& tgtMagSf
            ) const = 0;

            //- Normalise the weight. Can optionally subset addressing
            //- (e.g. for mapNearest)
            virtual void normaliseWeights
            (
                const bool verbose,
                AMIInterpolation<SourcePatch, TargetPatch>& inter
            ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeAMIMethod(AMIType)                                                 \
                                                                               \
    typedef AMIMethod<AMIType::sourcePatchType,AMIType::targetPatchType>       \
        AMIMethod##AMIType;                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(AMIMethod##AMIType, 0);                \
    defineTemplateRunTimeSelectionTable(AMIMethod##AMIType, components);


#define makeAMIMethodType(AMIType, Method)                                     \
                                                                               \
    typedef Method<AMIType::sourcePatchType,AMIType::targetPatchType>          \
        Method##AMIType;                                                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Method##AMIType, 0);                   \
                                                                               \
    AMIMethod<AMIType::sourcePatchType,AMIType::targetPatchType>::             \
        addcomponentsConstructorToTable<Method##AMIType>                       \
        add##Method##AMIType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AMIMethod.C"
    #include "AMIMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
