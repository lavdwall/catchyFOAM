/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::ReactingGSPhaseModel

Description
    Class which represents phases with volumetric reactions. Returns the
    reaction rate and heat.

SourceFiles
    ReactingGSPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingGSPhaseModel_H
#define ReactingGSPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ReactingGSPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class ReactionType>
class ReactingGSPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Reaction model
        autoPtr<ReactionType> reaction_;

        //- Gas exchange model enabled
        bool gasExchange_;
        bool chemistryT_;

        //- Gas phase name (only used when using gas exchange model)
        word gasPhaseName_;


public:

    // Constructors

        ReactingGSPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~ReactingGSPhaseModel();


    // Member Functions

        //- Correct the species mass fractions
        virtual void correctSpecies();

        //- Correct the reaction rates
        virtual void correctReactions();

        //- Return the fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Yi) const;

        //- Return the solid species reaction rates
        virtual tmp<fvScalarMatrix> Rs(volScalarField& Ysi) const;

        //- Return the gas species reaction rates
        virtual tmp<fvScalarMatrix> Rg(const volScalarField& Yi) const;

        //- Return heat release rate
        virtual tmp<volScalarField> Qdot() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingGSPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
