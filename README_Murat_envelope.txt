# README – Using `Murat_envelope.m`

## Overview

The function **`Murat_envelope.m`**, located in the `bin` directory, computes the envelope of seismic waveforms used during the MuRAT preprocessing workflow.

To correctly run this function, an additional parameter must be defined in the **MuRAT input file**.

## Required Input Parameter

In the **input configuration file** (e.g., `Murat_input.mlx` or similar), you must add the following line inside the **`Waveform Data`** section:


Murat.input.envelopeSmoothTime = 1;


## Parameter Description

* **`Murat.input.envelopeSmoothTime`**
  Defines the smoothing time (in seconds) applied when computing the waveform envelope.

* **Recommended value:**
  `1` second

This parameter is required by `Murat_envelope.m` to correctly smooth the envelope after it is computed from the seismic waveform.

## Where to Place the Parameter

Insert the line in the **Waveform Data** section of the input file, together with the other waveform-related parameters.

Example:


%% Waveform Data

Murat.input.envelopeSmoothTime = 1;


## Notes

* If this parameter is not defined, the envelope computation may fail or produce unexpected results.
* The value can be adjusted depending on the desired smoothing level, but **1 second is the standard value used in most MuRAT workflows**.

