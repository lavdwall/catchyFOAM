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
    Foam::SurfaceReaction

Description

    Read solid reactions of the type S1 = S2 + G1

SourceFiles
    SurfaceReactionI.H
    SurfaceReaction.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceReaction_H
#define SurfaceReaction_H

#include "speciesTable.H"
#include "specieCoeffs.H"
#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class ReactionThermo>
class SurfaceReaction;

template<class ReactionThermo>
inline Ostream& operator<<(Ostream&, const SurfaceReaction<ReactionThermo>&);

class objectRegistry;

/*---------------------------------------------------------------------------*\
                        Class SurfaceReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class SurfaceReaction
:
    public Reaction<ReactionThermo>
{

private:

    // Private data

        //- List of gas species present in reaction system
        speciesTable gases_;

        //- Number of gas species present in reaction system
        label nGases_;

        //- Gas specie index for the left-hand-side of the reaction
        List<specieCoeffs> glhs_;

        //- Gas specie index for the right-hand-side of the reaction
        List<specieCoeffs> grhs_;


    // Private Member Functions

        //- Return string representation of reaction
        string solidReactionStr(OStringStream&) const;

        //- Disallow default bitwise assignment
        void operator=(const SurfaceReaction&);


public:

    //- Runtime type information
    TypeName("SurfaceReaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            SurfaceReaction,
            dictionary,
            (
                const speciesTable& species,
                const HashPtrTable<ReactionThermo>& thermoDatabase,
                const dictionary& dict
            ),
            (species, thermoDatabase, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            SurfaceReaction,
            objectRegistry,
            (
                const speciesTable& species,
                const HashPtrTable<ReactionThermo>& thermoDatabase,
                const objectRegistry& ob,
                const dictionary& dict
            ),
            (species, thermoDatabase, ob, dict)
        );


    // Constructors

        //- Construct from components
        SurfaceReaction
        (
            const Reaction<ReactionThermo>& reaction,
            const speciesTable& gases,
            const List<specieCoeffs>& glhs,
            const List<specieCoeffs>& grhs
        );


        //- Construct as copy given new speciesTable
        SurfaceReaction
        (
            const SurfaceReaction<ReactionThermo>&,
            const speciesTable& gases
        );


        //- Construct from dictionary
        SurfaceReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );


    // Selectors
        static autoPtr<SurfaceReaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        static autoPtr<SurfaceReaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const objectRegistry& ob,
            const dictionary& dict
        );

        static autoPtr<SurfaceReaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const PtrList<ReactionThermo>& speciesThermo,
            const dictionary& dict
        );


    //- Destructor
    virtual ~SurfaceReaction()
    {}


    // Member Functions

        // Access

            //- Access to the gas components of the left hand side
            virtual const List<specieCoeffs>& grhs() const;

            //- Access to the gas components of the right hand side
            virtual const List<specieCoeffs>& glhs() const;

            //- Access to gas specie list
            virtual const speciesTable& gasSpecies() const;

        // Reaction rate coefficients

            //- Net reaction rate for individual species
            void omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            //- Net reaction rate
            scalar omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalar& pf,
                scalar& cf,
                label& lRef,
                scalar& pr,
                scalar& cr,
                label& rRef
            ) const;

        // Jacobian coefficients

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the species concentration
            void dwdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarSquareMatrix& J,
                scalarField& dcdt,
                scalar& omegaI,
                scalar& kfwd,
                scalar& kbwd,
                const bool reduced,
                const List<label>& c2s
            ) const;

        //- Coverage dependencies
        virtual inline const List<Tuple2<label, scalarField>>& covdep() const = 0;

        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionThermo>
        (
            Ostream&,
            const SurfaceReaction<ReactionThermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
