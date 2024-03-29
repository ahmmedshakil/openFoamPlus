/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2008-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "StandardWallInteraction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StandardWallInteraction<CloudType>::StandardWallInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    mesh_(cloud.mesh()),
    interactionType_
    (
        this->wordToInteractionType(this->coeffDict().getWord("type"))
    ),
    e_(0.0),
    mu_(0.0),
    nEscape_(mesh_.boundaryMesh().nNonProcessor()),
    massEscape_(nEscape_.size()),
    nStick_(nEscape_.size()),
    massStick_(nEscape_.size()),
    injIdToIndex_()
{
    const bool outputByInjectorId
        = this->coeffDict().lookupOrDefault("outputByInjectorId", false);

    switch (interactionType_)
    {
        case PatchInteractionModel<CloudType>::itOther:
        {
            const word interactionTypeName(this->coeffDict().getWord("type"));

            FatalErrorInFunction
                << "Unknown interaction result type "
                << interactionTypeName
                << ". Valid selections are:" << this->interactionTypeNames_
                << endl << exit(FatalError);

            break;
        }
        case PatchInteractionModel<CloudType>::itRebound:
        {
            e_ = this->coeffDict().lookupOrDefault("e", 1.0);
            mu_ = this->coeffDict().lookupOrDefault("mu", 0.0);
            break;
        }
        default:
        {}
    }

    // Determine the number of injectors and the injector mapping
    label nInjectors = 0;
    if (outputByInjectorId)
    {
        for (const auto& inj : cloud.injectors())
        {
            injIdToIndex_.insert(inj.injectorID(), nInjectors++);
        }
    }

    // The normal case, and safety if injector mapping was somehow null.
    if (injIdToIndex_.empty())
    {
        nInjectors = 1;
    }

    forAll(nEscape_, patchi)
    {
        nEscape_[patchi].setSize(nInjectors, Zero);
        massEscape_[patchi].setSize(nInjectors, Zero);
        nStick_[patchi].setSize(nInjectors, Zero);
        massStick_[patchi].setSize(nInjectors, Zero);
    }
}


template<class CloudType>
Foam::StandardWallInteraction<CloudType>::StandardWallInteraction
(
    const StandardWallInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    mesh_(pim.mesh_),
    interactionType_(pim.interactionType_),
    e_(pim.e_),
    mu_(pim.mu_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_),
    injIdToIndex_(pim.injIdToIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StandardWallInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    vector& U = p.U();

    if (isA<wallPolyPatch>(pp))
    {
        // Location for storing the stats.
        const label idx =
        (
            injIdToIndex_.size()
          ? injIdToIndex_.lookup(p.typeId(), 0)
          : 0
        );

        switch (interactionType_)
        {
            case PatchInteractionModel<CloudType>::itNone:
            {
                return false;
            }
            case PatchInteractionModel<CloudType>::itEscape:
            {
                keepParticle = false;
                p.active(false);
                U = Zero;

                const scalar dm = p.nParticle()*p.mass();

                nEscape_[pp.index()][idx]++;
                massEscape_[pp.index()][idx] += dm;
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                keepParticle = true;
                p.active(false);
                U = Zero;

                const scalar dm = p.nParticle()*p.mass();

                nStick_[pp.index()][idx]++;
                massStick_[pp.index()][idx] += dm;
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                p.active(true);

                vector nw;
                vector Up;

                this->owner().patchData(p, pp, nw, Up);

                // Calculate motion relative to patch velocity
                U -= Up;

                scalar Un = U & nw;
                vector Ut = U - Un*nw;

                if (Un > 0)
                {
                    U -= (1.0 + e_)*Un*nw;
                }

                U -= mu_*Ut;

                // Return velocity to global space
                U += Up;

                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type "
                    << this->interactionTypeToWord(interactionType_)
                    << "(" << interactionType_ << ")" << endl
                    << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


template<class CloudType>
void Foam::StandardWallInteraction<CloudType>::info(Ostream& os)
{
    PatchInteractionModel<CloudType>::info(os);

    labelListList npe0(nEscape_);
    this->getModelProperty("nEscape", npe0);

    scalarListList mpe0(massEscape_);
    this->getModelProperty("massEscape", mpe0);

    labelListList nps0(nStick_);
    this->getModelProperty("nStick", nps0);

    scalarListList mps0(massStick_);
    this->getModelProperty("massStick", mps0);

    // Accumulate current data
    labelListList npe(nEscape_);

    forAll(npe, i)
    {
        Pstream::listCombineGather(npe[i], plusEqOp<label>());
        npe[i] = npe[i] + npe0[i];
    }

    scalarListList mpe(massEscape_);
    forAll(mpe, i)
    {
        Pstream::listCombineGather(mpe[i], plusEqOp<scalar>());
        mpe[i] = mpe[i] + mpe0[i];
    }

    labelListList nps(nStick_);
    forAll(nps, i)
    {
        Pstream::listCombineGather(nps[i], plusEqOp<label>());
        nps[i] = nps[i] + nps0[i];
    }

    scalarListList mps(massStick_);
    forAll(nps, i)
    {
        Pstream::listCombineGather(mps[i], plusEqOp<scalar>());
        mps[i] = mps[i] + mps0[i];
    }

    if (injIdToIndex_.size())
    {
        // Since injIdToIndex_ is a one-to-one mapping (starting as zero),
        // can simply invert it.
        labelList indexToInjector(injIdToIndex_.size());
        forAllConstIters(injIdToIndex_, iter)
        {
            indexToInjector[iter.val()] = iter.key();
        }

        forAll(npe, i)
        {
            forAll(mpe[i], idx)
            {
                os  << "    Parcel fate: patch " <<  mesh_.boundary()[i].name()
                    << " (number, mass)" << nl
                    << "      - escape  (injector " << indexToInjector[idx]
                    << ")  = " << npe[i][idx]
                    << ", " << mpe[i][idx] << nl
                    << "      - stick   (injector " << indexToInjector[idx]
                    << ")  = " << nps[i][idx]
                    << ", " << mps[i][idx] << nl;
            }
        }
    }
    else
    {
        forAll(npe, i)
        {
            os  << "    Parcel fate: patch (number, mass) "
                << mesh_.boundary()[i].name() << nl
                << "      - escape                      = "
                << npe[i][0] << ", " << mpe[i][0] << nl
                << "      - stick                       = "
                << nps[i][0] << ", " << mps[i][0] << nl;

        }
    }

    if (this->writeTime())
    {
        this->setModelProperty("nEscape", npe);
        this->setModelProperty("massEscape", mpe);
        this->setModelProperty("nStick", nps);
        this->setModelProperty("massStick", mps);

        nEscape_ = Zero;
        massEscape_ = Zero;
        nStick_ = Zero;
        massStick_ = Zero;
    }
}


// ************************************************************************* //
