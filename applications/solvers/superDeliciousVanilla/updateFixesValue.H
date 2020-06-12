// If using timeVaryingMappedInletOutlet BCs, switch fixesValue to true (the
// boundary cannot be used to adjust the massflow to satisfy continuity) or 
// false (massflow is adjustable). Assume fixesValue should be true if over
// 50% of the flux is inward.
//
// Note: the actual continuity correction takes place in 
// src/finiteVolume/cfdTools/general/adjustPhi/adjustPhi.C

if(updateFixesValue)
{
    surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();

    forAll(bphi, patchi)
    {
        scalar massIn = 0.0;
        scalar massOut = 0.0;
        bool fixesValue;
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvsPatchScalarField& phip = bphi[patchi];

        if (isA<timeVaryingMappedInletOutletFvPatchField<vector> >(Up))
        {
            forAll(phip, i)
            {
                if (phip[i] < 0.0)
                {
                    massIn -= phip[i];
                }
                else
                {
                    massOut += phip[i];
                }
            }
            reduce(massIn, sumOp<scalar>());
            reduce(massOut, sumOp<scalar>());

            // verbose output
            //Info<< "[TVMIO] patch " << Up.patch().name()
            //    << " massIn= " << massIn
            //    << " massOut= " << massOut
            //    << endl;

            if(massIn > massOut)
            {
                fixesValue = true;
            }
            else
            {
                fixesValue = false;
            }

            if(fixesValue != Up.fixesValue())
            {
                Info<< "[TVMIO] switch " << Up.patch().name()
                    << " to fixesValue=" << fixesValue
                    << " at t=" << runTime.timeName()
                    << endl;

                /* Is this the best way to access the TVMIO patch? Reference:

                template<class TestType, class Type>
                inline bool isA(const Type& t)
                {
                    const Type* tPtr = &t;
                    return dynamic_cast<const TestType*>(tPtr);
                }

                */
                timeVaryingMappedInletOutletFvPatchField<vector>* tvmioPatch =
                    dynamic_cast<timeVaryingMappedInletOutletFvPatchField<vector>*>
                    (
                        &U.boundaryFieldRef()[patchi]
                    );

                tvmioPatch->setFixesValue(fixesValue);
            }
        }
    }
}
