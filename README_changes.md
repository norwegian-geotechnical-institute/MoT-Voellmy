# MoT-Voellmy – Change History

<br/>
This document collects the changes between versions, starting from v2024-09-10.

### Changes from v2026-04-20 to v2026-05-04

* Fixed bug in the entrainment model `AvaFrame`: Entrainment (or comminution energy is to be read from file, not from number `Erosion coefficient` (in function read_command_file).
* Validity check of specific comminution energy (entire raster) in function read_init corrected.
* Documentation adjusted to new version
* Additional RCF for Flateyri avalanche, showcasing the `AvaFrame` entrainment formula
* Version number in Makefile updated


### Changes from v2025-05-20 to v2026-04-20

* New RCF (Run Control File) format 2026-04-20:
    - Line 3 shows the version date 2026-04-20.
    - On line 59, the variable `Momentum threshold` (in kg m/s) is replaced with `Minimum momentum fraction` (in %).
* `MoT-Voellmy v2026-04-20` is compatible with SCF versions `2020-06-23`, `2021-10-25`, `2024-09-10`, and `2026-04-20`.
* Older versions since 2021-10-25 allowed to stop the simulation if the integrated "quantity of movement" *J* = ∫ *ρ h* |**u**| d*A* falls below a predefined value *J*<sub>min</sub> set by the user. In v2026-04-20, the user specifies a *fraction* of the maximum quantity of movement, *J*<sub>max</sub>, attained in the simulation as the stopping criterion.
* The erosion energy for the AvaFrame entrainment model is set to 5 kJ/kg if the user specified a value outside the range recommended by BFW (0.2–10 kJ/kg).
* Slightly reduced output to `stdout` (typically the screen or terminal).


### Changes from v2024-09-10 to v2025-05-20

* Restored a factor |**u**|² in the calculation of the centrifugal force, which had been erroneously omitted in v2024-09-10. Thanks are due to Hervé Vicari for spotting the error.


