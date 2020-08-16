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
    Foam::SurfaceReactionProxy

Description
    Proxy version of SolidReaction which provides dummy implementations of the
    abstract virtual functions.

    Used for read-construction and format conversion.

SourceFiles
    SurfaceReactionProxyI.H
    SurfaceReactionProxy.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceReactionProxy_H
#define SurfaceReactionProxy_H

#include "SurfaceReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SurfaceReactionProxy Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class SurfaceReactionProxy
:
    public SurfaceReaction<ReactionThermo>
{

public:

    // Constructors

        //- Construct from components
        SurfaceReactionProxy
        (
            const Reaction<ReactionThermo>& reaction,
            const speciesTable& gases,
            const List<specieCoeffs>& glhs,
            const List<specieCoeffs>& grhs
        );


        //- Construct as copy given new speciesTable
        SurfaceReactionProxy
        (
            const SurfaceReactionProxy<ReactionThermo>&,
            const speciesTable& gases
        );


        //- Construct from dictionary
        SurfaceReactionProxy
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ReactionThermo>> clone() const;

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const;


    //- Destructor
    virtual ~SurfaceReactionProxy()
    {}


    // Member Functions

        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Reverse rate constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;


        // Jacobian coefficients

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Temperature derivative of reverse rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar dkfdT,
                const scalar kr
            ) const;

            //- Third-body efficiencies (beta = 1-alpha)
            //  non-empty only for third-body reactions
            //  with enhanced molecularity (alpha != 1)
            virtual const List<Tuple2<label, scalar>>& beta() const;

            //- Coverage dependencies
            virtual inline const List<Tuple2<label, scalarField>>& covdep() const;

            //- Species concentration derivative of the pressure dependent term
            virtual void dcidc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dcidc
            ) const;

            //- Temperature derivative of the pressure dependent term
            virtual scalar dcidT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceReactionProxy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
