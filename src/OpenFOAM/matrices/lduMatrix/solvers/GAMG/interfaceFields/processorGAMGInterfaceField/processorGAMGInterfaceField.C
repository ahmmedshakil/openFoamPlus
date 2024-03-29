/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2019 OpenCFD Ltd.
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

#include "processorGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        processorGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        processorGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorGAMGInterfaceField::processorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    procInterface_(refCast<const processorGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const processorLduInterfaceField& p =
        refCast<const processorLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::processorGAMGInterfaceField::processorGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    procInterface_(refCast<const processorGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorGAMGInterfaceField::~processorGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorGAMGInterfaceField::initInterfaceMatrixUpdate
(
    solveScalarField&,
    const bool,
    const solveScalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    procInterface_.interfaceInternalField(psiInternal, scalarSendBuf_);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        scalarReceiveBuf_.setSize(scalarSendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        IPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<char*>(scalarReceiveBuf_.begin()),
            scalarReceiveBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        OPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<const char*>(scalarSendBuf_.begin()),
            scalarSendBuf_.byteSize(),
            procInterface_.tag(),
            comm()
        );
    }
    else
    {
        procInterface_.compressedSend(commsType, scalarSendBuf_);
    }

    const_cast<processorGAMGInterfaceField&>(*this).updatedMatrix() = false;
}


void Foam::processorGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const solveScalarField&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (updatedMatrix())
    {
        return;
    }

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }
        // Recv finished so assume sending finished as well.
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;

        // Consume straight from scalarReceiveBuf_

        // Transform according to the transformation tensor
        transformCoupleField(scalarReceiveBuf_, cmpt);

        // Multiply the field by coefficients and add into the result
        addToInternalField(result, !add, coeffs, scalarReceiveBuf_);
    }
    else
    {
        solveScalarField pnf
        (
            procInterface_.compressedReceive<solveScalar>
            (
                commsType,
                coeffs.size()
            )
        );
        transformCoupleField(pnf, cmpt);

        addToInternalField(result, !add, coeffs, pnf);
    }

    const_cast<processorGAMGInterfaceField&>(*this).updatedMatrix() = true;
}


// ************************************************************************* //
