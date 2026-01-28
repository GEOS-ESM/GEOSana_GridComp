# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

28Jan2026:
- bring back changes from PR #209

28Jan2026:
- bug fix: add ABI tlapmean back (had lost them when backing out PR#209)

12Jan2026:
- bug fixes in interfacing GSI w/ libSP - these are zero-diff fixes

09Jan2026:
- Back out changes from PR 209 (IR adjustments) since this is 
  aimed at FPP and the update must be zero-diff.
- optional write out in setupw

07Jan2026:
- minor revisions of CG tolerance (consistency between pcgsoi and bicg)

04Sep2025:

- add knobs for high-resolution radiosondes (zero-diff
  when new obsclass not used).

16May2025:
 - revise handling of OMPS-LP for 5.42.x FPP/FP

29Apr2025:
- add settings to allow proper comparison with JEDI.

02Apr2025:
- Revisions to original FSI implemetation:
   o apply it every outer iteration
   o leave original omf alone

12Mar2025:
- add FSI knob (off by default)
- set ta2tb as true by default (in main config files)
- set number of crtm surfaces to 20 as in presently in GEOS-JEDI

12Dec2024:
- updates GeoVals falling the Oct 2024 JEDI variable name change code sprint.
- bug fix in setupbend to allow running old flavor of pressure/integration option.
- add operationally-named class for OMPS-LP from NPP has been added to gsi.rc.tmpl and gsi_sens.rc.tmpl

------------------------------------------------------
- add ability to convert antenna to brigtness temperature
- bugfix: update calculation of div/vor
- merge of 5.30.3 with 5.40.0 and 5.29.5-p7
- add offline RO bufr handling program
- add CrIS-FSR N21
- revert obs errors for CrIS-Npp and N20 to what they
  were in x0049.
- add fix to read_bufrtovs to handle ta2tb=.true. when
  there are multiple versions of SpcCoeff.bin file for 
  single instrument/platform
- add ability to handle nvege_type=20 - no longer doing
  fishy interpolation.
- update to use CRTM v2.4.1-jedi-1
- Revise sources of OMPS-LP; fix handling of multiple 
  sources of OMPS-LP in source code.
- minor bug fix in gsi_sens.rc.tmpl for CrIS/N21 dsis


### Fixed

### Removed

### Deprecated

