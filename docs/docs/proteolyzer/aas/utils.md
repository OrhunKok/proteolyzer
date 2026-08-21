---
sidebar_label: utils
title: proteolyzer.aas.utils
---

Reference tables and small helpers shared by the AAS pipeline stages.

## Path

## np

## pd

## reference

## write\_frame

#### column\_mapping

```python
def column_mapping(df: pd.DataFrame, cols2keep: list) -> pd.DataFrame
```

Standardized column selection + renaming based on &#x27;cols2keep&#x27;.

Matching is case-insensitive, and the output column order follows
``cols2keep`` rather than the order of the input file.

#### aa\_subs\_ref

```python
def aa_subs_ref() -> dict
```

Mass delta of every single amino acid substitution, keyed by origin residue.

``{&quot;A&quot;: {&quot;A to G&quot;: -14.0157, ...}, ...}``

#### gen\_mod\_dict

```python
def gen_mod_dict() -> dict
```

Known modifications per residue as ``[full_name, position, mono_mass]`` rows.

Excludes the ``AA substitution`` class, which is what :func:`aa_subs_ref`
returns. This reference exists to decide whether a mass shift has an
explanation *other than* mistranslation, so leaving the substitutions in it
would make every candidate substitution match itself as a modification.

#### saap\_alt\_output

```python
def saap_alt_output(dp_df: pd.DataFrame, sample_name: str,
                    output_dir: Path) -> None
```

Split detected dependent peptides into the ALT and SAAP sets.

A SAAP carries an amino acid substitution that no known modification
explains (those go to ALT) and that has no homologue in any translated
frame. BASE is the unmodified peptide a SAAP is measured against.

#### calculate\_aa\_substitution\_matrix

```python
def calculate_aa_substitution_matrix(
        processed_amino_acids_df: pd.DataFrame) -> pd.DataFrame
```

Calculates the pairwise mass difference matrix (Row AA mass - Column AA mass).

