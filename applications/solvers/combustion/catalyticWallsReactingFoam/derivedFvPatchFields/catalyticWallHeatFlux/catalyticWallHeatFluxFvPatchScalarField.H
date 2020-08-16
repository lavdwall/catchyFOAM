/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::catalyticWallHeatFluxFvPatchScalarField

Description


Usage


See also
    Foam::fixedGradientFvPatchField

SourceFiles
    catalyticWallHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef catalyticWallHeatFluxFvPatchScalarField_H
#define catalyticWallHeatFluxFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class catalyticWallHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class catalyticWallHeatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
protected:

    //- Catalyst name
    const word catalystName_;

public:

    //- Runtime type information
    TypeName("catalyticWallHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        catalyticWallHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        catalyticWallHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given catalyticWallHeatFluxFvPatchScalarField
        // onto a new patch
        catalyticWallHeatFluxFvPatchScalarField
        (
            const catalyticWallHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        catalyticWallHeatFluxFvPatchScalarField
        (
            const catalyticWallHeatFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new catalyticWallHeatFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        catalyticWallHeatFluxFvPatchScalarField
        (
            const catalyticWallHeatFluxFvPatchScalarField&,
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
                new catalyticWallHeatFluxFvPatchScalarField(*this, iF)
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
