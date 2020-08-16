/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::psiReactionThermo

Description
    Foam::psiReactionThermo

SourceFiles
    psiReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiReactionThermo_H
#define psiReactionThermo_H

#include "psiThermo.H"
#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class psiReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class psiReactionThermo
:
    public psiThermo
{

    //- List of active reacting cells
    List<bool> reactingCells_;

public:

    //- Runtime type information
    TypeName("psiReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        psiReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiReactionThermo();


    // Member Functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;

        //- Set reacting status of cell, celli
        virtual void setCellReacting(const label celli, const bool active);

        //- Get reacting status of cell, celli
        bool cellReacting(const label celli) const {return reactingCells_[celli];}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
