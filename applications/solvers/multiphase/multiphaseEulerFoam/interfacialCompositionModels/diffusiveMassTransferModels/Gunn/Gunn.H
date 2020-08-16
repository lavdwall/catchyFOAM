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
    Foam::diffusiveMassTransferModels::Gunn

Description

SourceFiles
    Gunn.C

\*---------------------------------------------------------------------------*/

#ifndef Gunn_H
#define Gunn_H

#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace diffusiveMassTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class Gunn Declaration
\*---------------------------------------------------------------------------*/

class Gunn
:
    public diffusiveMassTransferModel
{
    // Private Data

        //- Lewis number
        const dimensionedScalar Le_;


public:

    //- Runtime type information
    TypeName("Gunn");


    // Constructors

        //- Construct from components
        Gunn
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Gunn();


    // Member Functions

        //- The implicit mass transfer coefficient
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusiveMassTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
