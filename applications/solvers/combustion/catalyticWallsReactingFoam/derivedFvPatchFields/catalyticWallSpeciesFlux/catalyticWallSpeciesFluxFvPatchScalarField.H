/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::catalyticWallSpeciesFluxFvPatchScalarField

Description


Usage


See also
    Foam::fixedGradientFvPatchField

SourceFiles
    catalyticWallSpeciesFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef catalyticWallSpeciesFluxFvPatchScalarField_H
#define catalyticWallSpeciesFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class catalyticWallSpeciesFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class catalyticWallSpeciesFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
protected:

    //- Catalyst name
    const word catalystName_;

public:

    //- Runtime type information
    TypeName("catalyticWallSpeciesFlux");


    // Constructors

        //- Construct from patch and internal field
        catalyticWallSpeciesFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        catalyticWallSpeciesFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given catalyticWallSpeciesFluxFvPatchScalarField
        // onto a new patch
        catalyticWallSpeciesFluxFvPatchScalarField
        (
            const catalyticWallSpeciesFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        catalyticWallSpeciesFluxFvPatchScalarField
        (
            const catalyticWallSpeciesFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new catalyticWallSpeciesFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        catalyticWallSpeciesFluxFvPatchScalarField
        (
            const catalyticWallSpeciesFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new catalyticWallSpeciesFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
