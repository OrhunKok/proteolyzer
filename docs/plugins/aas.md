# AAS: substitution discovery

The pipeline used for discovery of amino acid substitutions, from the
[Slavov lab](https://decode.slavovlab.net/). Every stage reads the same
parameter file — see `examples/aas/params.yaml`.

## Nomenclature

| term | meaning |
|---|---|
| **SAAP** | a peptide carrying an amino acid substitution |
| **BASE** | the unmodified peptide a SAAP is measured against |
| **ALT** | a mass shift with an alternative explanation: a known modification |

Outputs live under `SAAP/` and `ALT/`, and the quantification columns are
`SAAP.Sum`, `BASE.Sum`, `SAAP.Plex.<n>`, `BASE.Plex.<n>`. The mass tolerance
for matching a known modification is `ALT ppm`.

## Two phases

There is a manual database search in the middle: detection writes a FASTA that
has to be searched against the raw files before validation has anything to
read.

```python
import proteolyzer.aas as aas

pipeline = aas.Pipeline("params.yaml")

pipeline.run_detection()    # preprocess -> translate -> detect
#   ... search the raw files against <output>/<sample>_validation.fasta ...
pipeline.run_validation()   # preprocess -> validate -> quantify

pipeline.status()           # what has run, and what can run now
```

`Pipeline` exists for the ordering: the preprocessor runs again in phase two to
convert the validation searches, the six-frame translation is skipped when its
frames are already on disk, and phase two refuses to start before the searches
exist.

The stages can still be driven individually:

```python
aas.Preprocessor.MaxQuant(params).run()
aas.FrameTranslator(params).run()
aas.Detection(params).run()
aas.Validation(params).run()
aas.Quantification(params).run()
```

Stages exchange frames as parquet under the output folder, and each run appends
its resolved parameters, timestamp and package version to
`<output folder>/provenance.jsonl`.

## Reading a run back

Without knowing the stage-internal file names:

```python
results = aas.Results.from_params(params)   # or aas.Results("out/")

results.samples                # what is in there
results.summary()              # rows per sample per step; NA where a step did
                               # not run, so reading down a column shows where
                               # it stopped
results.combined("quantified") # every sample in one frame, with a Sample column
results.provenance()           # what produced it
```

Steps are named for their result rather than the file that holds it:
`candidates`, `alt`, `filtered`, `fasta_entries`, `validated`, `evidence`,
`quantified`.

!!! note "MaxQuant only"

    `Utils.Workflow` accepts MaxQuant. The stages are structured so a second
    search engine slots in alongside, but nothing else is implemented.
