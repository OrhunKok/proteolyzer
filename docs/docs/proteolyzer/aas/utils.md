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

#### ptm\_mtp\_output

```python
def ptm_mtp_output(dp_df: pd.DataFrame, sample_name: str,
                   output_dir: Path) -> None
```

Split detected dependent peptides into PTM and MTP sets and pickle both.

MTP are peptides with a potential AA substitution that cannot be explained
by a PTM and that have no homologous sequence in any translated frame.

#### calculate\_aa\_substitution\_matrix

```python
def calculate_aa_substitution_matrix(
        processed_amino_acids_df: pd.DataFrame) -> pd.DataFrame
```

Calculates the pairwise mass difference matrix (Row AA mass - Column AA mass).

