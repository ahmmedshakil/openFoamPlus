/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Modified code Copyright (C) 2013-2019 FOSS GP
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

#include "adjointWallVelocityLowReFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointWallVelocityLowReFvPatchVectorField::
adjointWallVelocityLowReFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointBoundaryCondition(p, iF, "Ua")
{}


Foam::adjointWallVelocityLowReFvPatchVectorField::
adjointWallVelocityLowReFvPatchVectorField
(
    const adjointWallVelocityLowReFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    adjointBoundaryCondition(p, iF, "Ua")
{}


Foam::adjointWallVelocityLowReFvPatchVectorField::
adjointWallVelocityLowReFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointBoundaryCondition(p, iF, "Ua")
{
    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
}


Foam::adjointWallVelocityLowReFvPatchVectorField::
adjointWallVelocityLowReFvPatchVectorField
(
    const adjointWallVelocityLowReFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    adjointBoundaryCondition(pivpvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointWallVelocityLowReFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Objective function contribution
    tmp<vectorField> tsource = boundaryContrPtr_->normalVelocitySource();
    vectorField& source = tsource.ref();

    operator==(-source);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointWallVelocityLowReFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
    os.writeEntry("solverName", adjointSolverName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointWallVelocityLowReFvPatchVectorField
    );
}


// ************************************************************************* //
