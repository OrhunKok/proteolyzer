# Reference data and UniMod

## Constants

Masses, the genetic code and protease rules live in one place, as immutable
mappings read on first use:

```python
from proteolyzer import reference

reference.amino_acid_masses()["A"]          # 71.03711...
reference.protease("Trypsin").allowed_counts
reference.CODON_TABLE["AUG"]                # "M"
reference.three_letter_to_one()["Ala"]      # "A"
```

`reference` is imported with the package and holds no settings — only facts
about the domain. Anything configurable lives in
[`core.formats`][proteolyzer.core.formats] instead.

## Querying UniMod

The two bundled tables are exported from UniMod. For anything they cannot
answer, query the database itself:

```python
from proteolyzer import unimod

unimod.tables()                 # built on first use
unimod.table("modifications")

unimod.query(
    """
    SELECT m.full_name, s.one_letter, m.mono_mass
    FROM specificity AS s
    JOIN modifications AS m ON s.mod_key = m.record_id
    WHERE m.mono_mass BETWEEN ? AND ?
    """,
    (15.9, 16.1),
)
```

The database is 5 MB against 221 KB of CSVs, so it is not shipped: it is built
once into `$PROTEOLYZER_CACHE_DIR` (else `$XDG_CACHE_HOME`, else `~/.cache`)
and reused.

**Querying uses only the standard library.** Building downloads the UniMod XML
and needs `pip install '.[unimod]'` — which is why the query API and the
builder are separate.

## Regenerating the bundled tables

A maintainer step, after a new UniMod release:

```bash
python -m proteolyzer.unimod build
python -m proteolyzer.unimod export \
    --mods-output src/proteolyzer/resources/unimod_modifications.csv \
    --aa-output   src/proteolyzer/resources/unimod_amino_acids.csv
```
