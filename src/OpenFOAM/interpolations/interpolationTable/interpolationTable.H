/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2004-2011 OpenCFD Ltd.
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::interpolationTable

Description
    An interpolation/look-up table of scalar vs \<Type\> values.
    The reference scalar values must be monotonically increasing.

    The handling of out-of-bounds values depends on the current setting
    of \c outOfBounds.

    If \c repeat is chosen for the out-of-bounds handling, the final time
    value is treated as being equivalent to time=0 for the following periods.


    The construct from dictionary reads a filename from a dictionary and
    has an optional readerType. Default is to read OpenFOAM format. The only
    other format is csv (comma separated values):

    Read csv format:
    \verbatim
        readerType      csv;
        file            "<constant>/p0vsTime.csv";
        hasHeaderLine   true;   // skip first line
        timeColumn      0;      // time is in column 0
        valueColumns    (1);    // value starts in column 1
    \endverbatim


Note
    - Accessing an empty list results in an error.
    - Accessing a list with a single element always returns the same value.

SourceFiles
    interpolationTable.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationTable_H
#define interpolationTable_H

#include "List.H"
#include "Tuple2.H"

#include "tableBounds.H"
#include "tableReader.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interpolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationTable
:
    public List<Tuple2<scalar, Type>>
{
    // Private data

        //- Handling for out-of-bound values
        bounds::repeatableBounding bounding_;

        //- File name
        fileName fileName_;

        //- The actual reader
        autoPtr<tableReader<Type>> reader_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();


public:

    // Constructors

        //- Construct null
        interpolationTable();

        //- Construct from components
        interpolationTable
        (
            const List<Tuple2<scalar, Type>>& values,
            const bounds::repeatableBounding bounding,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        interpolationTable(const fileName& fName);

        //- Construct by reading file name and outOfBounds from dictionary
        //  and read the table from that file.
        //  This is a specialised constructor used by patchFields
        interpolationTable(const dictionary& dict);

        //- Construct copy
        interpolationTable(const interpolationTable& interpTable);


    // Member Functions

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;

        //- Write
        void write(Ostream& os) const;

        //- Return the rate of change at the interpolation location
        //  for the give value
        Type rateOfChange(const scalar) const;


    // Member Operators

        //- Return an element of constant Tuple2<scalar, Type>
        const Tuple2<scalar, Type>& operator[](const label) const;

        //- Return an interpolated value
        Type operator()(const scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
