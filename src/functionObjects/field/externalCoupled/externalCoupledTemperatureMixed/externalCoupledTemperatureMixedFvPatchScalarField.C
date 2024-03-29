/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "externalCoupledTemperatureMixedFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "Enum.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::externalCoupledTemperatureMixedFvPatchScalarField::
    outputTemperatureType
>
Foam::externalCoupledTemperatureMixedFvPatchScalarField::outputTemperatureNames
({
    { outputTemperatureType::FLUID, "fluid" },
    { outputTemperatureType::WALL, "wall" },
});


const Foam::Enum
<
    Foam::externalCoupledTemperatureMixedFvPatchScalarField::
    refTemperatureType
>
Foam::externalCoupledTemperatureMixedFvPatchScalarField::refTemperatureNames
({
    { refTemperatureType::CELL, "cell" },
    { refTemperatureType::USER, "user" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeHeader
(
    Ostream& os
) const
{
    if (outTempType_ == outputTemperatureType::WALL)
    {
        os  << "# Values: area Twall qDot htc" << endl;
    }
    else
    {
        os  << "# Values: area Tfluid qDot htc" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(p, iF),
    outTempType_(outputTemperatureType::WALL),
    refTempType_(refTemperatureType::CELL),
    Tref_(Zero)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    externalCoupledMixedFvPatchField<scalar>(ptf, p, iF, mapper),
    outTempType_(ptf.outTempType_),
    refTempType_(ptf.refTempType_),
    Tref_(ptf.Tref_)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    //externalCoupledMixedFvPatchField<scalar>(p, iF, dict)
    externalCoupledMixedFvPatchField<scalar>(p, iF),
    outTempType_(outputTemperatureType::WALL),
    refTempType_
    (
        refTemperatureNames.lookupOrDefault
        (
            "htcRefTemperature",
            dict,
            refTemperatureType::CELL
        )
    ),
    Tref_(Zero)
{
    if (dict.found("outputTemperature"))
    {
        outTempType_ = outputTemperatureNames.get("outputTemperature", dict);
    }
    else
    {
        WarningInFunction
            << "outputTemperature not specified "
            << flatOutput(outputTemperatureNames) << nl
            << "using 'wall' as compatibility default" << nl
            << endl;
    }

    if (refTempType_ == refTemperatureType::USER)
    {
        Tref_ = dict.get<scalar>("Tref");
    }

    if (dict.found("refValue"))
    {
        // Initialise same way as mixed
        this->refValue() = scalarField("refValue", dict, p.size());
        this->refGrad() = scalarField("refGradient", dict, p.size());
        this->valueFraction() = scalarField("valueFraction", dict, p.size());

        evaluate();
    }
    else
    {
        // For convenience: initialise as fixedValue with either read value
        // or extrapolated value
        if (dict.found("value"))
        {
            fvPatchField<scalar>::operator=
            (
                scalarField("value", dict, p.size())
            );
        }
        else
        {
            fvPatchField<scalar>::operator=(this->patchInternalField());
        }

        // Initialise as a fixed value
        this->refValue() = *this;
        this->refGrad() = Zero;
        this->valueFraction() = 1.0;
    }
}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf
)
:
    externalCoupledMixedFvPatchField<scalar>(ecmpf),
    outTempType_(ecmpf.outTempType_),
    refTempType_(ecmpf.refTempType_),
    Tref_(ecmpf.Tref_)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(ecmpf, iF),
    outTempType_(ecmpf.outTempType_),
    refTempType_(ecmpf.refTempType_),
    Tref_(ecmpf.Tref_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeData
(
    Ostream& os
) const
{
    const label patchi = patch().index();

    // Heat flux [W/m2]
    scalarField qDot(this->patch().size(), Zero);

    typedef compressible::turbulenceModel cmpTurbModelType;

    static word turbName
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    static word thermoName("thermophysicalProperties");

    if (db().foundObject<cmpTurbModelType>(turbName))
    {
        const cmpTurbModelType& turbModel =
            db().lookupObject<cmpTurbModelType>(turbName);

        const basicThermo& thermo = turbModel.transport();

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchi];

        qDot = turbModel.alphaEff(patchi)*hep.snGrad();
    }
    else if (db().foundObject<basicThermo>(thermoName))
    {
        const basicThermo& thermo = db().lookupObject<basicThermo>(thermoName);

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchi];

        qDot = thermo.alpha().boundaryField()[patchi]*hep.snGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Condition requires either compressible turbulence and/or "
            << "thermo model to be available" << exit(FatalError);
    }


    // Wall temperature [K]
    const scalarField& Twall = *this;

    // Fluid temperature [K]
    tmp<scalarField> tfluid;

    if (refTempType_ == refTemperatureType::USER)
    {
        // User-specified reference temperature
        tfluid = tmp<scalarField>::New(size(), Tref_);
    }
    else
    {
        // Near wall cell temperature
        tfluid = patchInternalField();
    }

    const scalarField Tfluid(tfluid);


    // Heat transfer coefficient [W/m2/K]
    const scalarField htc(qDot/(max(Twall - Tfluid), 1e-3));

    const Field<scalar>& magSf = this->patch().magSf();

    const UList<scalar>& Tout =
    (
        outTempType_ == outputTemperatureType::FLUID
      ? Tfluid
      : Twall
    );

    forAll(patch(), facei)
    {
        os  << magSf[facei] << token::SPACE
            << Tout[facei] << token::SPACE
            << qDot[facei] << token::SPACE
            << htc[facei] << nl;
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::readData
(
    Istream& is
)
{
    // Assume generic input stream so we can do line-based format and skip
    // unused columns
    ISstream& iss = dynamic_cast<ISstream&>(is);

    string line;

    forAll(*this, facei)
    {
        iss.getLine(line);
        IStringStream lineStr(line);

        lineStr
            >> this->refValue()[facei]
            >> this->refGrad()[facei]
            >> this->valueFraction()[facei];
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    externalCoupledMixedFvPatchField::write(os);
    os.writeEntry
    (
        "outputTemperature",
        outputTemperatureNames[outTempType_]
    );
    os.writeEntry
    (
        "htcRefTemperature",
        refTemperatureNames[refTempType_]
    );

    if (refTempType_ == refTemperatureType::USER)
    {
        os.writeEntry("Tref", Tref_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalCoupledTemperatureMixedFvPatchScalarField
    );
}


// ************************************************************************* //
