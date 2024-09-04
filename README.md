# AffineMotions

[![Build Status](https://github.com/olivierverdier/AffineMotions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/olivierverdier/AffineMotions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/olivierverdier/AffineMotions.jl/graph/badge.svg?token=aTe2GSxvIw)](https://codecov.io/gh/olivierverdier/AffineMotions.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierverdier.github.io/AffineMotions.jl/)

Given a group action $G \subset \mathrm{Diff}(M)$
a function $φ \colon M \to G$ is *affine* if the quantity
```math
φχ(χ\cdot x) - χ φ(x) χ^{-1}
```
does not depend on $x$.

This package defines various operations for those functions.
