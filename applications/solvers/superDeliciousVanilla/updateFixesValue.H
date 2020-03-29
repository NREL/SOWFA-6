// If using timeVaryingMappedInletOutlet BCs, switch fixesValue to true (the
// boundary cannot be used to adjust the massflow to satisfy continuity) or 
// false (massflow is adjustable). Assume fixesValue should be true if over
// 50% of the flux is inward.
//
// Note: the actual continuity correction takes place in 
// src/finiteVolume/cfdTools/general/adjustPhi/adjustPhi.C
if(updateFixesValue)
{
    Info<< "Updating fixesValue" << endl;
}
