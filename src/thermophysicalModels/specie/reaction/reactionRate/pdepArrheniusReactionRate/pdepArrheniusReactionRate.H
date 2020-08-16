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
    Foam::pdepArrheniusReactionRate

Description
    General pressure dependence using logarithmic interpolation of Arrhenius
    expressions.

SourceFiles
    pdepArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef pdepArrheniusReactionRate_H
#define pdepArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class pdepArrheniusReactionRate;

Ostream& operator<<(Ostream&, const pdepArrheniusReactionRate&);


/*---------------------------------------------------------------------------*\
                    Class pdepArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class pdepArrheniusReactionRate
{
    // Private data

        scalarList p_;
        scalarList A_;
        scalarList beta_;
        scalarList Ta_;


public:

    // Constructors

        //- Construct from components
        inline pdepArrheniusReactionRate
        (
            const scalarList p,
            const scalarList A,
            const scalarList beta,
            const scalarList Ta
        );

        //- Construct from dictionary
        inline pdepArrheniusReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "pdepArrhenius";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label
        ) const;

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        inline const List<Tuple2<label, scalar>>& beta() const;

        //- Species concentration derivative of the pressure dependent term
        //  By default this value is 1 as it multiplies the third-body term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
        //  By default this value is 0 since ddT of molecularity is approx.0
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const pdepArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdepArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
