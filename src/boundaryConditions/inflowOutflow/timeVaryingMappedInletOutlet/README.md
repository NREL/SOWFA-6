# Changes from v2.4.x

## Changelog

- 2021-06-20 - fixed offset and setAverage bug
- 2020-06-03 - initial implementation in OpenFOAM-6


## Details

- Constructor from patch, internal field, and dictionary: If a value is provided in the BC dict,
  then apply operator== instead of operator= (change in OpenFOAM-6). Otherwise, call update the
  fixed values and then apply operator== instead of calling evaluate(). Evaluate() calls:
  ```
  void Foam::mixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
  {
      if (!this->updated())
      {
          this->updateCoeffs();
      }

      Field<Type>::operator=
      (
          valueFraction_*refValue_
        +
          (1.0 - valueFraction_)*
          (
              this->patchInternalField()
            + refGrad_/this->patch().deltaCoeffs()
          )
      );

      fvPatchField<Type>::evaluate();
  }
  ```
  However, this is a non-starter because the call to updateCoeffs() will fail prior to the first
  time step because phi has not been calculated. Instead, we call updateFixedValue()--the original
  part of timeVaryingMappedFixedValue--followed by operator== to set the value to the provided
  reference values from constant/boundaryData. _This may have been avoided in the past by always
  specifying the value keyword._

- Separate updateCoeffs into separate, self-explanatory `updateFixedValue()` and
  `updateInletOutlet()` functions.

- `valueFraction` is calculated as 1 - pos0(phi), where pos0 is 1 if >= 0 and 0 otherwise, instead
  of 1 - pos(phi) (OpenFOAM-6 change).

- `evaluate()` is now overridden to provide statistics about valueFraction.

- `operator=` is no longer overridden; this was a redefinition of the more general calculation in
  evaluate(), but without the "refGradient" term. In our case, refGradient is always set to zero 
  and the inherited operator= and operator== should be consistent. 
  
- It is important to remember that this `TVMIO` class is a child of the `mixed` boundary condition 
  class, not the `TVMFV` class.  This code is basically a copy of `TVMFV` with the necessary 
  changes and inheritance of `mixed` instead of `fixedValue`.  As such, it is important to remember
  that the code should not be applying the boundary data values of quantities to the actual `value`
  of the patch, but rather to `refVal()`.  Once `refVal` is set, `evaluate()` in the `mixed` boundary
  condition class will compute what the patch values should be based on whether it is inflow or
  outflow.  Therefore, when the code reads in boundary data, that data should be used to set `refValue()`.


## Minor changes affecting usage

- Added `format` keyword (either "ascii" or "binary") to allow input points and sample files to
  be in uncompressed/compressed ascii/binary format.
- `AverageIOField` has been replaced in OpenFOAM-6 by `AverageField`, which is no longer derived 
  from the `regIOobject` class. The "FoamFile" header dictionary should therefore be omitted.
- Added optional keywords "dataDir", "points", and "sample" in OpenFOAM-6.
