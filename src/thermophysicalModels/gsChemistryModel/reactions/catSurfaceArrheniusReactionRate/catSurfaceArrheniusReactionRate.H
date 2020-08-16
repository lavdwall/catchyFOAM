/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::catSurfaceArrheniusReactionRate

Description
    Arrhenius reaction rate for solids

SourceFiles
    catSurfaceArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef catSurfaceArrheniusReactionRate_H
#define catSurfaceArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class catSurfaceArrheniusReactionRate;

Ostream& operator<<(Ostream&, const catSurfaceArrheniusReactionRate&);


/*---------------------------------------------------------------------------*\
                       Class catSurfaceArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class catSurfaceArrheniusReactionRate
{
    // Private data

        scalar A_;
        scalar beta_;
        scalar Ta_;


public:

    // Constructors

        //- Construct from components
        inline catSurfaceArrheniusReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta
        );


        //- Construct from dictionary
        inline catSurfaceArrheniusReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    virtual ~catSurfaceArrheniusReactionRate()
    {}


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "SurfaceArrhenius";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Third body efficiencies
        inline const List<Tuple2<label, scalar>>& beta() const;

        //- Coverage dependencies
        inline const List<Tuple2<label, scalarField>>& covdep() const;

        //- Species concentration derivative of the coverage dependent term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the coverage dependent term
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;


         //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const catSurfaceArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "catSurfaceArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
