---
sidebar_label: detection
title: proteolyzer.aas.detection
---

Detection of candidate amino acid substitutions and PTMs.

:class:`Detection` is a thin dispatcher: it reads the ``Utils.Workflow``
parameter and delegates to the matching subclass (currently only
:class:`MaxQuant`).

## pickle

## re

## shutil

## ahocorasick

## np

## pd

## reference

## read\_frame

## write\_frame

## utils

## Stage

## Config

#### CONFIG

## Detection Objects

```python
class Detection(Stage)
```

Dispatches to the search-engine specific detection workflow.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### \_initialize\_workflow

```python
def _initialize_workflow()
```

Instantiate the subclass named by the configured workflow.

#### run

```python
def run()
```

Run the detection workflow.

#### process\_sample

```python
def process_sample(sample)
```

Delegate sample processing to the detection workflow.

## MaxQuant Objects

```python
class MaxQuant(Detection)
```

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### run

```python
def run()
```

Run the MaxQuant-specific detection workflow over every sample.

#### process\_sample

```python
def process_sample(sample)
```

Process each sample with MaxQuant-specific logic.

#### \_validate\_saap

```python
def _validate_saap(sample, evidence_path)
```

Validate SAAPs and write the filtered set.

#### find\_potential\_aas

```python
def find_potential_aas(dp_df, n_term_peps, c_term_peps)
```

Identify potential AA substitutions and other modifications.

#### \_prepare\_peptide\_data

```python
def _prepare_peptide_data(dp_df, n_term_peps, c_term_peps)
```

Prepare the peptide data for further processing.

#### \_apply\_modifications

```python
def _apply_modifications(dp_df)
```

Assign substitutions, then alternative explanations, to peptides.

#### refine\_localization\_probabilities

```python
def refine_localization_probabilities(modified_seq)
```

#### pep\_cterm

```python
def pep_cterm(modified_sequence)
```

#### pep\_nterm

```python
def pep_nterm(modified_sequence)
```

#### \_find\_homologous\_peptides

```python
def _find_homologous_peptides(dp_df)
```

Find homologous peptides based on translated frames.

#### get\_aa\_subs

```python
def get_aa_subs(row)
```

Get the amino acid substitutions for a given row.

#### get\_saap\_sequence

```python
def get_saap_sequence(row)
```

Generate the SAAP sequence implied by the substitution.

#### find\_alt

```python
def find_alt(row)
```

Find an alternative explanation (a known modification) for a peptide.

#### aho\_corasick\_search

```python
def aho_corasick_search(translated_genome, patterns)
```

Search for patterns in the translated genome using Aho-Corasick.

#### aho\_corasick\_output\_organize

```python
def aho_corasick_output_organize(result)
```

Organize the output from the Aho-Corasick search.

#### find\_homologous\_peptide

```python
def find_homologous_peptide(sample_df)
```

Find homologous peptides using the Aho-Corasick algorithm.

#### evidence\_ppm

```python
def evidence_ppm(evidence)
```

Calculate the average and standard deviation of mass error (ppm).

#### posterior\_aasub\_prob

```python
def posterior_aasub_prob(row, mean, std)
```

Calculate posterior probability for amino acid substitution.

#### q\_val\_calc

```python
def q_val_calc(saap)
```

Calculate q-values for peptides based on posterior probabilities.

#### gen\_metrics

```python
def gen_metrics(saap)
```

Generate evaluation metrics.

#### saap\_filter

```python
def saap_filter(saap, metric_df)
```

Filter SAAPs based on F-score.

#### validate\_saap

```python
def validate_saap(evidence, saap)
```

Validate SAAPs using posterior probability and q-value filtering.

#### write\_fasta

```python
def write_fasta(filtered_saap, sample_name)
```

Write the filtered SAAPs to a FASTA file.

