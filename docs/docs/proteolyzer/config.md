---
sidebar_label: config
title: proteolyzer.config
---

Domain reference data shared by every proteolyzer subpackage.

Amino acid masses and three-to-one letter codes are read once from the bundled
UniMod export in ``proteolyzer/resources``; see :mod:`proteolyzer.unimod` for
how those CSVs are regenerated.

## resources

## io

## dataclass

## Traversable

## ClassVar

## Final

## pd

#### \_RESOURCE\_MODULE

#### \_AA\_FILE

#### \_get\_resource\_stream

```python
def _get_resource_stream(filename: str) -> io.BytesIO
```

## AminoAcids Objects

```python
@dataclass(frozen=True)
class AminoAcids()
```

Container for amino acid properties.

#### MASS

#### CODE\_TO\_SYMBOL

## Codons Objects

```python
@dataclass(frozen=True)
class Codons()
```

Container for genetic code mappings, organized by species.

#### Standard

fmt: skip

## Protease Objects

```python
@dataclass(frozen=True)
class Protease()
```

Defines cleavage rules and allowed residue counts for common proteases.
Each protease is a nested subclass with:
    - CLEAVAGE_SITES: residues where cleavage occurs
    - ALLOWED_COUNTS: maximum allowed counts of specific residues in peptides

CLEAVAGE_SITES are tuples: they are class-level shared state, so callers
must build a new sequence rather than extending them in place.

## Trypsin Objects

```python
@dataclass(frozen=True)
class Trypsin()
```

#### CLEAVAGE\_SITES

#### ALLOWED\_COUNTS

## LysC Objects

```python
@dataclass(frozen=True)
class LysC()
```

#### CLEAVAGE\_SITES

#### ALLOWED\_COUNTS

## ArgC Objects

```python
@dataclass(frozen=True)
class ArgC()
```

#### CLEAVAGE\_SITES

#### ALLOWED\_COUNTS

