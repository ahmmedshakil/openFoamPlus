/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::mapNearestAMI

Description
    Nearest-mapping Arbitrary Mesh Interface (AMI) method

SourceFiles
    mapNearestAMI.C

\*---------------------------------------------------------------------------*/

#ifndef mapNearestAMI_H
#define mapNearestAMI_H

#include "AMIMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mapNearestAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class mapNearestAMI
:
    public AMIMethod<SourcePatch, TargetPatch>
{

private:

    // Private Member Functions

        //- No copy construct
        mapNearestAMI(const mapNearestAMI&) = delete;

        //- No copy assignment
        void operator=(const mapNearestAMI&) = delete;

        // Marching front

            //- Find nearest target face for source face srcFacei
            void findNearestFace
            (
                const SourcePatch& srcPatch,
                const TargetPatch& tgtPatch,
                const label& srcFacei,
                label& tgtFacei
            ) const;

            //- Determine next source-target face pair
            void setNextNearestFaces
            (
                boolList& mapFlag,
                label& startSeedI,
                label& srcFacei,
                label& tgtFacei
            ) const;

            //- Find mapped source face
            label findMappedSrcFace
            (
                const label tgtFacei,
                const List<DynamicList<label>>& tgtToSrc
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFacei,
                const label tgtFacei
            ) const;


public:

    //- Runtime type information
    TypeName("mapNearestAMI");


    // Constructors

        //- Construct from components
        mapNearestAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );


    //- Destructor
    virtual ~mapNearestAMI();


    // Member Functions

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
            );

            //- Set the face areas for parallel runs
            virtual void setMagSf
            (
                const TargetPatch& tgtPatch,
                const mapDistribute& map,
                scalarList& srcMagSf,
                scalarList& tgtMagSf
            ) const;

            //- Normalise the weight. Can optionally subset addressing
            //- (e.g. for mapNearest)
            virtual void normaliseWeights
            (
                const bool verbose,
                AMIInterpolation<SourcePatch, TargetPatch>& inter
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mapNearestAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
