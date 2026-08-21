---
sidebar_label: aas
title: proteolyzer.aas
---

Amino acid substitution (AAS) discovery pipeline.

Implements the pipeline used for discovery of amino acid substitutions and
PTMs (https://decode.slavovlab.net/). The stages are meant to be run in order,
each reading the same parameter file:

    Preprocessor -&gt; FrameTranslator -&gt; Detection -&gt; Validation -&gt; Quantification

Modules
    base: shared stage plumbing (parameter loading, queue, metadata)
    params: parameter file schema and loader
    results: read back what a run produced
    preprocessing: search-engine output to parquet conversion
    translation: six-frame translation of a genome FASTA
    detection: SAAP and ALT assignment
    validation: fragment-level validation of candidates
    quantification: substitution ratio quantification

## NullQueue

## Stage

## Detection

## MaxQuant

## ParamsSchema

## load\_params

## Preprocessor

## Quantification

## ARTEFACTS

## Results

## FrameTranslator

## Validation

#### \_\_all\_\_

