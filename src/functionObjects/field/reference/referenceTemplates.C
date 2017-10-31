/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "interpolation.H"

template<class Type>
bool Foam::functionObjects::reference::calcType()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* vfPtr = lookupObjectPtr<VolFieldType>(fieldName_);

    if (vfPtr)
    {
        const VolFieldType& vf = *vfPtr;

        dimensioned<Type> offset
        (
            dimensioned<Type>::lookupOrDefault
            (
                "offset",
                localDict_,
                vf.dimensions(),
                Zero
            )
        );

        dimensioned<Type> cellValue("value", vf.dimensions(), Zero);

        if (positionIsSet_)
        {
            cellValue.value() = -pTraits<Type>::one*GREAT;

            if (celli_ != -1)
            {
                autoPtr<interpolation<Type>> interpolator
                (
                    interpolation<Type>::New(interpolationScheme_, vf)
                );

                cellValue.value() =
                    interpolator().interpolate(position_, celli_, -1);
            }

            reduce(cellValue.value(), maxOp<Type>());

            Log << "    sampled value: " << cellValue.value() << endl;
        }

        return store
        (
            resultName_,
            scale_*(vf - cellValue + offset)
        );
    }

    return false;
}


// ************************************************************************* //