---
sidebar_label: reference
title: proteolyzer.reference
---

Domain reference data: one source of truth for the constants.

Everything here is fixed knowledge about proteins and proteases rather than a
setting, so it is exposed as module-level immutable mappings and cached
loaders instead of a configuration object. Subpackages import what they need
directly; they do not re-expose it through their own settings objects, which
only added a layer to look through.

The tables come from the UniMod export in ``proteolyzer/resources``, read on
first use rather than at import, so importing the package costs nothing. See
``tools/unimod`` for how those files are regenerated.

## functools

## resources

## io

## Mapping

## dataclass

## MappingProxyType

## Final

## pd

#### \_RESOURCE\_MODULE

#### AMINO\_ACID\_FILE

#### MODIFICATION\_FILE

#### resource\_stream

```python
def resource_stream(filename: str) -> io.BytesIO
```

The bundled resource `filename`, as a stream.

#### amino\_acids

```python
@functools.cache
def amino_acids() -> pd.DataFrame
```

The amino acid table, including Xle (``J``) for isobaric Ile/Leu.

#### modifications

```python
@functools.cache
def modifications() -> pd.DataFrame
```

Every approved UniMod modification, one row per residue it applies to.

#### amino\_acid\_masses

```python
@functools.cache
def amino_acid_masses() -> Mapping[str, float]
```

Monoisotopic residue mass by one-letter code.

#### three\_letter\_to\_one

```python
@functools.cache
def three_letter_to_one() -> Mapping[str, str]
```

One-letter code by three-letter code, e.g. ``{&quot;Ala&quot;: &quot;A&quot;}``.

#### \_frozen

```python
def _frozen(mapping: dict) -> Mapping
```

A read-only view, so a shared table cannot be mutated by a caller.

#### CODON\_TABLE

The standard genetic code, keyed by RNA codon. ``*`` is a stop codon.

## Protease Objects

```python
@dataclass(frozen=True)
class Protease()
```

Cleavage rules for one protease.

Attributes
----------
cleavage_sites
    Residues the protease cleaves after.
allowed_counts
    How many times each of those residues may appear in a fully cleaved
    peptide. A tuple and a read-only mapping, because these are shared:
    extending them in place used to corrupt every later call.

#### name

#### cleavage\_sites

#### allowed\_counts

#### PROTEASES

Proteases the processing pipeline can flag missed cleavages for.

#### protease

```python
def protease(name: str) -> Protease
```

The named protease, or a ValueError listing the ones there are.

