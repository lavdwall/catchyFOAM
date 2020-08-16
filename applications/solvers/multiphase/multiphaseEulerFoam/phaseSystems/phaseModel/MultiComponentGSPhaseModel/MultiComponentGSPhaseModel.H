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
    Foam::MultiComponentGSPhaseModel

Description
    Class which represents a phase with multiple species. Returns the species'
    mass fractions, and their governing equations.

SourceFiles
    MultiComponentGSPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MultiComponentGSPhaseModel_H
#define MultiComponentGSPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MultiComponentGSPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MultiComponentGSPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- Inert species index
        label inertIndex_;

        //- Inert solid species index
        label inertSIndex_;

        //- Pointer list to active species
        UPtrList<volScalarField> YActive_;

        //- Pointer list to active solid species
        UPtrList<volScalarField> YsActive_;


public:

    // Constructors

        MultiComponentGSPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~MultiComponentGSPhaseModel();


    // Member Functions

        //- Correct the species fractions
        virtual void correctSpecies();

        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi);

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef();

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef();

            //- Return the solid species fraction equation
            virtual tmp<fvScalarMatrix> YsiEqn(volScalarField& Ysi);

            //- Return the solid species mass fractions
            virtual const PtrList<volScalarField>& Ys() const;

            //- Return a solid species mass fraction by name
            virtual const volScalarField& Ys(const word& name) const;

            //- Access the solid species mass fractions
            virtual PtrList<volScalarField>& YsRef();

            //- Return the active solid species mass fractions
            virtual const UPtrList<volScalarField>& YsActive() const;

            //- Access the active solid species mass fractions
            virtual UPtrList<volScalarField>& YsActiveRef();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MultiComponentGSPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
