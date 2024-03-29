/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Released 2009-2011 OpenCFD Ltd.
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Modified code Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::functionObjects::fieldAverage

Group
    grpFieldFunctionObjects

Description
    Calculates average quantities for a user-specified selection of volumetric
    and surface fields.
    With the %subRegion option, also supports fields on functionObject
    surface output (eg, sampledSurfaces).


    Fields are entered as a list of sub-dictionaries, which indicate the type of
    averages to perform, and can be updated during the calculation.  The current
    options include:
    - \c mean: arithmetic mean
        \f[
            \overline{x} = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N x_i
        \f]
    - \c prime2Mean: prime-squared mean
        \f[
            \overline{x'}^2 = \frac{1}{N}\displaystyle\sum\limits_{i=0}^N
            (x_i - \overline{x})^2
        \f]
    - \c base: average over 'time', or 'iteration' (\f$N\f$ in the above)
    - \c window: optional averaging window, specified in 'base' units

    Average field names are constructed by concatenating the base field with
    the averaging type, e.g. when averaging field 'U', the resultant fields
    are:
    - arithmetic mean field, \c UMean
    - prime-squared field, \c UPrime2Mean

    Information regarding the number of averaging steps, and total averaging
    time are written on a per-field basis to the
    \c "<functionObject name>Properties" dictionary, located in \<time\>/uniform

    When restarting form a previous calculation, the averaging is continuous or
    may be restarted using the \c restartOnRestart option.

    The averaging process may be restarted after each calculation output time
    using the \c restartOnOutput option or restarted periodically using the \c
    periodicRestart option and setting \c restartPeriod to the required
    averaging period.

Usage
    \verbatim
    fieldAverage1
    {
        type                fieldAverage;
        libs                ("libfieldFunctionObjects.so");

        writeControl        writeTime;

        restartOnRestart    false;
        restartOnOutput     false;
        periodicRestart     false;
        restartPeriod       0.002;

        fields
        (
            U
            {
                mean            on;
                prime2Mean      on;
                base            time;
                windowType      exact;
                window          10.0;
                windowName      w1;
            }
            p
            {
                mean            on;
                prime2Mean      on;
                base            time;
            }
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property        | Description                       | Required | Default
        type            | Type name: fieldAverage           | yes |
        restartOnRestart| Restart the averaging on restart  | no  | no
        restartOnOutput | Restart the averaging on output   | no  | no
        periodicRestart | Periodically restart the averaging | no  | no
        restartPeriod   | Periodic restart period            | conditional |
        restartTime     | One-shot reset of the averaging    | no | great
        fields          | list of fields and averaging options | yes |
        subRegion       | Name for alternative objectRegistry  | no  | ""
    \endtable


Note
    To employ the \c prime2Mean option, the \c mean option must be selecetd.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    fieldAverage.C
    fieldAverageTemplates.C
    fieldAverageItem.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldAverage_H
#define functionObjects_fieldAverage_H

#include "fvMeshFunctionObject.H"
#include "FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward declaration of classes
class fieldAverageItem;

/*---------------------------------------------------------------------------*\
                        Class fieldAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldAverage
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Time at last call, prevents repeated averaging
        label prevTimeIndex_;

        //- Initialised flag
        bool initialised_;

        //- Restart the averaging process on restart
        Switch restartOnRestart_;

        //- Restart the averaging process on output
        Switch restartOnOutput_;

        //- Periodically restart the averaging process
        Switch periodicRestart_;

        //- Restart period
        scalar restartPeriod_;

        //- Specific restart time
        scalar restartTime_;

        //- List of field average items, describing what averages to be
        //  calculated and output
        List<fieldAverageItem> faItems_;

        // Counters

            //- Iteration steps counter
            List<label> totalIter_;

            //- Total time counter
            List<scalar> totalTime_;

            //- Index for periodic restart
            label periodIndex_;


    // Protected Member Functions

        // Initialisation routines

            //- Reset lists (clear existing values) and initialize averaging.
            //  Check requested field averages are valid, populate field lists
            void initialize();

            //- Restart averaging for restartOnOutput
            void restart();

            //- Add mean average field to database
            template<class Type>
            void addMeanFieldType(fieldAverageItem& item);

            //- Add mean average field to database
            template<class Type>
            void addMeanField(fieldAverageItem& item);

            //- Add prime-squared average field to database
            template<class Type1, class Type2>
            void addPrime2MeanFieldType(fieldAverageItem& item);

            //- Add prime-squared average field to database
            template<class Type1, class Type2>
            void addPrime2MeanField(fieldAverageItem& item);


        // Calculation functions

            //- Main calculation routine
            virtual void calcAverages();

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFields() const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFields() const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2MeanType(const fieldAverageItem& item) const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2Mean() const;

            template<class Type>
            void storeWindowFieldType(fieldAverageItem& item);

            template<class Type>
            void storeWindowFields();

            template<class Type>
            void restoreWindowFieldsType(const fieldAverageItem& item);

            template<class Type>
            void restoreWindowFields(const fieldAverageItem& item);

        // I-O

            //- Write averages
            virtual void writeAverages() const;

            //- Write fields
            template<class Type>
            void writeFieldType(const word& fieldName) const;

            //- Write fields
            template<class Type>
            void writeFields() const;

            //- Write averaging properties - steps and time
            void writeAveragingProperties();

            //- Read averaging properties - steps and time
            void readAveragingProperties();


        //- No copy construct
        fieldAverage(const fieldAverage&) = delete;

        //- No copy assignment
        void operator=(const fieldAverage&) = delete;


public:

    //- Runtime type information
    TypeName("fieldAverage");


    // Constructors

        //- Construct from Time and dictionary
        fieldAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~fieldAverage() = default;


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
