# QPhiX Wrapper
A wrapper around the [QPhiX library](https://github.com/JeffersonLab/qphix) for lattice QCD developed by Jefferson Lab to simplify calls to the QPhiX solvers

## Description
QPhiX is a highly optimized library for performing propagator inversions for lattice QCD calculations on CPU architectures with vectorized SIMD instructions.  It was developed by Jefferson Lab as part of the USQCD collaboration, and major lattice QCD software packages (e.g. chroma) have the ability to link against QPhiX's optimized inverters to improve performance.

The most optimized inverters in the QPhiX library assume even-odd preconditioning, which reduces the condition number of the matrix being inverted (and thus the computational effort required) by a factor of 2 or more.  Typically, this preconditioning is performed by the calling software package (such as chroma).

The motivation of this package is to provide a standalone, encapsulated solver for an unpreconditioned fermion that performs the preconditioning, the inversion, and the reconstruction of both even and odd sites.  In particular, the goal was a single function that would be given a fermion $\psi$ and return $D^{-1} \psi$.  In addition, this function should have the following properties:
- callable from both C and C++ (and in particular hides the C++ objects from C code that would call it)
- stores all internal variables needed by both the inverter and the preconditioner to avoid expensive regeneration of inverse clover terms, etc.
- codebase should be minimal and fast to install

The codebase here is based heavily on the clover inversion test included in QPhiX, and most of the codebase is copied from the testing scripts.  Since the purpose of this wrapper is not testing but production, it exposes functions so that fermion and gauge fields can be passed in and out.

Many thanks to Balint Joo for very helpful discussions explaining the mechanics of the QPhiX preconditioner that were essential in understanding and replicating the behavior of a traditional chroma/QPhiX setup.

## Instructions

## Branches

In addition to the main branch, two others are worth noting:
- `0vbb`: This branch was used for a calculation of neutrinoless double-beta decay and is called by our codebase used for [that project](https://github.com/agrebe/0vbb).  It has been frozen in time to avoid breaking linkages.
- `qc`: This is the primary development branch.  It is designed to be used for a separate lattice QCD codebase called QC that is still in development and has not yet been publicly released.  This branch may have additional features but is also less stable (and may contain breaking changes).
