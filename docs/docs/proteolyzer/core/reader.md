---
sidebar_label: reader
title: proteolyzer.core.reader
---

The entry point: turn a path into a :class:`Report`.

## Data

## Report

#### read

```python
def read(source, **kwargs) -> Report
```

Read search-engine output into a :class:`Report`.

``source`` is a path, a string path, or a file-like object. Remaining
keyword arguments configure :class:`Data`, e.g. ``load_all_columns=True``
or ``extra_cols_to_load=[&quot;Q.Value&quot;]``.

    report = proteolyzer.read(&quot;report.parquet&quot;)
    processed = report.process()

