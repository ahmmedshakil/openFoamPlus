/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2008-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "ReactingMultiphaseParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingMultiphaseParcel<ParcelType>::propertyList_ =
    Foam::ReactingMultiphaseParcel<ParcelType>::propertyList();

template<class ParcelType>
Foam::string Foam::ReactingMultiphaseParcel<ParcelType>::propertyTypes_ =
    Foam::ReactingMultiphaseParcel<ParcelType>::propertyTypes();

template<class ParcelType>
const std::size_t Foam::ReactingMultiphaseParcel<ParcelType>::sizeofFields
(
    0
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0)
{
    if (readFields)
    {
        DynamicList<scalar> Yg;
        DynamicList<scalar> Yl;
        DynamicList<scalar> Ys;

        is >> Yg >> Yl >> Ys;

        YGas_.transfer(Yg);
        YLiquid_.transfer(Yl);
        YSolid_.transfer(Ys);

        // scale the mass fractions
        const scalarField& YMix = this->Y_;
        YGas_ /= YMix[GAS] + ROOTVSMALL;
        YLiquid_ /= YMix[LIQ] + ROOTVSMALL;
        YSolid_ /= YMix[SLD] + ROOTVSMALL;
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    bool valid = c.size();

    ParcelType::readFields(c, compModel);

    // Get names and sizes for each Y...
    const label idGas = compModel.idGas();
    const wordList& gasNames = compModel.componentNames(idGas);
    const label idLiquid = compModel.idLiquid();
    const wordList& liquidNames = compModel.componentNames(idLiquid);
    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);
    const wordList& stateLabels = compModel.stateLabels();

    // Set storage for each Y... for each parcel
    for (ReactingMultiphaseParcel<ParcelType>& p : c)
    {
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
    }

    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas
        (
            c.fieldIOobject
            (
                "Y" + gasNames[j] + stateLabels[idGas],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (ReactingMultiphaseParcel<ParcelType>& p : c)
        {
            p.YGas_[j] = YGas[i]/(p.Y()[GAS] + ROOTVSMALL);
            ++i;
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject
            (
                "Y" + liquidNames[j] + stateLabels[idLiquid],
                 IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (ReactingMultiphaseParcel<ParcelType>& p : c)
        {
            p.YLiquid_[j] = YLiquid[i]/(p.Y()[LIQ] + ROOTVSMALL);
            ++i;
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j] + stateLabels[idSolid],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (ReactingMultiphaseParcel<ParcelType>& p : c)
        {
            p.YSolid_[j] = YSolid[i]/(p.Y()[SLD] + ROOTVSMALL);
            ++i;
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();

    // Write the composition fractions
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject
                (
                    "Y" + gasNames[j] + stateLabels[idGas],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            for (const ReactingMultiphaseParcel<ParcelType>& p0 : c)
            {
                YGas[i] = p0.YGas()[j]*p0.Y()[GAS];
                ++i;
            }

            YGas.write(np > 0);
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject
                (
                    "Y" + liquidNames[j] + stateLabels[idLiquid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            for (const ReactingMultiphaseParcel<ParcelType>& p0 : c)
            {
                YLiquid[i] = p0.YLiquid()[j]*p0.Y()[LIQ];
                ++i;
            }

            YLiquid.write(np > 0);
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject
                (
                    "Y" + solidNames[j] + stateLabels[idSolid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            for (const ReactingMultiphaseParcel<ParcelType>& p0 : c)
            {
                YSolid[i] = p0.YSolid()[j]*p0.Y()[SLD];
                ++i;
            }

            YSolid.write(np > 0);
        }
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    const CompositionType& compModel,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    label np = c.size();

    // Write the composition fractions
    if (np > 0)
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            const word fieldName = "Y" + gasNames[j] + stateLabels[idGas];
            IOField<scalar>& YGas
            (
                cloud::createIOField<scalar>(fieldName, np, obr)
            );

            label i = 0;
            for (const ReactingMultiphaseParcel<ParcelType>& p0 : c)
            {
                YGas[i] = p0.YGas()[j]*p0.Y()[GAS];
                ++i;
            }
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            const word fieldName = "Y" + liquidNames[j] + stateLabels[idLiquid];
            IOField<scalar>& YLiquid
            (
                cloud::createIOField<scalar>(fieldName, np, obr)
            );

            label i = 0;
            for (const ReactingMultiphaseParcel<ParcelType>& p0 : c)
            {
                YLiquid[i] = p0.YLiquid()[j]*p0.Y()[LIQ];
                ++i;
            }
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            const word fieldName = "Y" + solidNames[j] + stateLabels[idSolid];
            IOField<scalar>& YSolid
            (
                cloud::createIOField<scalar>(fieldName, np, obr)
            );

            label i = 0;
            for (const ReactingMultiphaseParcel<ParcelType>& p0 : c)
            {
                YSolid[i] = p0.YSolid()[j]*p0.Y()[SLD];
                ++i;
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingMultiphaseParcel<ParcelType>& p
)
{
    scalarField YGasLoc(p.YGas()*p.Y()[0]);
    scalarField YLiquidLoc(p.YLiquid()*p.Y()[1]);
    scalarField YSolidLoc(p.YSolid()*p.Y()[2]);
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << YGasLoc
            << token::SPACE << YLiquidLoc
            << token::SPACE << YSolidLoc;
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os  << YGasLoc << YLiquidLoc << YSolidLoc;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
