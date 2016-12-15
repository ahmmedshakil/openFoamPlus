/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

Description
    Containers for holding ccm solution and field listings.

\*---------------------------------------------------------------------------*/
#ifndef ccmSolutionTable_H
#define ccmSolutionTable_H

#include "SLList.H"
#include "stringListOps.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

/*---------------------------------------------------------------------------*\
                      Class ccm::namesList Declaration
\*---------------------------------------------------------------------------*/

//- A linked-list that is searchable by the 'name()' of the items
template<class T>
class namesList
:
    public SLList<T>
{
public:
    typedef typename SLList<T>::const_iterator const_iterator;
    typedef typename SLList<T>::iterator iterator;

    // Constructors

        //- Null construct
        namesList()
        {}


    // Access

        //- Return true if a list element has a name that matches key
        bool found(const word& key) const
        {
            for
            (
                const_iterator iter = SLList<T>::begin();
                iter != SLList<T>::end();
                ++iter
            )
            {
                if (iter().name() == key)
                {
                    return true;
                }
            }

            return false;
        }


        //- Find a list element has a name matching key
        iterator find(const word& key)
        {
            for
            (
                iterator iter = SLList<T>::begin();
                iter != SLList<T>::end();
                ++iter
            )
            {
                if (iter().name() == key)
                {
                    return iter;
                }
            }

            return SLList<T>::end();
        }


        //- Return a list of names matching whiteList and not matching blackList
        List<word> findNames
        (
            const UList<wordRe>& whiteLst,
            const UList<wordRe>& blackLst = UList<wordRe>()
        ) const
        {
            List<word> matched(SLList<T>::size());

            label matchI = 0;
            for
            (
                const_iterator iter = SLList<T>::begin();
                iter != SLList<T>::end();
                ++iter
            )
            {
                const word& name = iter().name();

                if
                (
                    findStrings(whiteLst, name)
                 && !findStrings(blackLst, name)
                )
                {
                    matched[matchI++] = name;
                }
            }
            matched.setSize(matchI);

            return matched;
        }

};


/*---------------------------------------------------------------------------*\
                      Class ccm::fieldEntry Declaration
\*---------------------------------------------------------------------------*/

//- A ccm field entry with short name, name, maxId and type
//  shortName => ( fullName, maxId, type );
class fieldEntry
{
    // Private Data

        //- The field name (PROSTAR short name)
        word name_;

        //- The full field name
        string fullName_;

        //- The field units
        string units_;

        //- The max cell id for the field
        label maxCellId_;

        //- The max face id for the field
        label maxFaceId_;

public:

    // Constructors

        //- Construct from components with optional units
        fieldEntry
        (
            const word& shortName,
            const string& fullName,
            const char* units = nullptr
        )
        :
            name_(shortName),
            fullName_(fullName),
            units_(),
            maxCellId_(0),
            maxFaceId_(0)
        {
            if (units && *units)
            {
                units_ = units;
            }
        }


    // Access

        //- The field name (PROSTAR short name)
        const word& name() const
        {
            return name_;
        }

        //- The full field name
        const string& fullName() const
        {
            return fullName_;
        }

        //- The field units
        const string& units() const
        {
            return units_;
        }

        //- The max cell id for the field
        label maxCellId() const
        {
            return maxCellId_;
        }

        //- The max face id for the field
        label maxFaceId() const
        {
            return maxFaceId_;
        }


    // Edit

        //- Set the field units
        void units(const char* units)
        {
            if (units && *units)
            {
                units_ = units;
            }
        }

        //- Set the field units
        void units(const std::string& units)
        {
            if (!units.empty())
            {
                units_ = units;
            }
        }

        //- Set the max cell Id for the field
        void maxCellId(const int newMax)
        {
            if (maxCellId_ < newMax)
            {
                maxCellId_ = newMax;
            }
        }


        //- Set the max face Id for the field
        void maxFaceId(const int newMax)
        {
            if (maxFaceId_ < newMax)
            {
                maxFaceId_ = newMax;
            }
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const fieldEntry& entry
        )
        {
            os  << entry.name_ << " => " << entry.fullName_
                << " [" << entry.units_.c_str()
                << "] maxCell: " << entry.maxCellId_
                << " maxFace: " << entry.maxFaceId_;

            return os;
        }

};

/*---------------------------------------------------------------------------*\
                    Class ccm::solutionEntry Declaration
\*---------------------------------------------------------------------------*/

//- A ccm solution entry with name, iteration and time
//  stateName => ( iteration, time );
class solutionEntry
{
    // Private Data

        //- The solution name
        word name_;

        //- The solution iteration/timestep
        label iter_;

        //- The solution time (sec)
        scalar time_;

public:

    // Constructors

        //- Construct from components
        solutionEntry
        (
            const word& name,
            const label& iteration,
            const scalar& timeValue = 0
        )
        :
            name_(name),
            iter_(iteration),
            time_(timeValue)
        {}

    // Access

        //- The solution name
        const word& name() const
        {
            return name_;
        }

        //- The solution iteration/timestep
        label iteration() const
        {
            return iter_;
        }

        //- The solution time (sec)
        scalar timeValue() const
        {
            return time_;
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const solutionEntry& entry
        )
        {
            os  << entry.name_ << " =>"
                << " iter: " << entry.iter_
                << " time: " << entry.time_;

            return os;
        }

};


/*---------------------------------------------------------------------------*\
                    Class ccm::solutionTable Declaration
\*---------------------------------------------------------------------------*/

// Typedef: ccm::solutionTable
// A list of all the available solutions
typedef namesList<solutionEntry> solutionTable;


/*---------------------------------------------------------------------------*\
                      Class ccm::fieldTable Declaration
\*---------------------------------------------------------------------------*/

//- A list of the available fields
class fieldTable
:
    public namesList<fieldEntry>
{
public:

    // Constructor

        //- Null construct
        fieldTable()
        :
            namesList<fieldEntry>()
        {}

    // Access

        //- The maximum cell Id referenced in the list
        label maxCellId() const
        {
            label maxId = 0;
            forAllConstIter(namesList<fieldEntry>, *this, iter)
            {
                label currMax = (iter()).maxCellId();

                if (maxId < currMax)
                {
                    maxId = currMax;
                }
            }

            return maxId;
        }


        //- The maximum face Id referenced in the list
        label maxFaceId() const
        {
            label maxId = 0;

            forAllConstIter(namesList<fieldEntry>, *this, iter)
            {
                label currMax = (iter()).maxFaceId();

                if (maxId < currMax)
                {
                    maxId = currMax;
                }
            }

            return maxId;
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const fieldTable& tbl
        )
        {
            os  << static_cast<const namesList<fieldEntry>& >(tbl)
                << nl
                << "maxCell: "  << tbl.maxCellId()
                << " maxFace: " << tbl.maxFaceId();

            return os;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //