# Chemistry configuration

The `tenx_chemistry` setting resolves the read geometry, barcode whitelist,
barcode translation, HAM chemistry, and UMI length used by guide extraction.
Definitions are stored in [`config/chemistry_spec.yaml`](../config/chemistry_spec.yaml).

```yaml
tenx_chemistry: "3v3"
```

Supported values are `3v3`, `3v4`, `3LT`, `multiome`, `5v1`, `5v2`, `5v3`,
and `custom`.

## 3′ direct capture

These chemistries use bead-borne cs1/cs2 RT primers and dual-oligo barcodes.
Guide and GEX libraries may use different barcode representations, so automatic
translation is required.

| Setting | 10x kit | R1 layout | UMI | Whitelist | Translation | Status |
|---|---|---|---:|---|:---:|:---:|
| `3v3` | 3′ v3/v3.1 | 16 bp CB + 12 bp UMI | 12 bp | `3M-feb-2018` | Yes | Validated |
| `3v4` | 3′ v4/GEM-X | 16 bp CB + 12 bp UMI | 12 bp | `3M-3pgex-may-2023` | Yes | Pending real-data validation |
| `3LT` | 3′ LT | 16 bp CB + 12 bp UMI | 12 bp | `3M-feb-2018` | Yes | Pending real-data validation |
| `multiome` | Multiome GEX | 16 bp CB + 12 bp UMI | 12 bp | `3M-feb-2018` | Yes | Pending real-data validation |

`3v3` and `3LT` share the same library structure. `3v4` uses a distinct
whitelist and translation table. The `3M-3pgex-may-2023` whitelist has no public
mirror and must be supplied from Cell Ranger 8.0.1 or later.

## 5′ direct capture

These chemistries use a soluble scaffold RT primer and a barcoded template
switch oligo. Guide and GEX libraries use the same barcode representation, so
translation is not performed.

| Setting | 10x kit | R1 layout | UMI | Whitelist | Translation | Status |
|---|---|---|---:|---|:---:|:---:|
| `5v1` | 5′ v1.0 | 16 bp CB + 10 bp UMI | 10 bp | `737K-aug-2016` | No | Validated |
| `5v2` | 5′ v2 | 16 bp CB + 10 bp UMI | 10 bp | `737K-aug-2016` | No | Pending real-data validation |
| `5v3` | 5′ v3/GEM-X | 16 bp CB + 12 bp UMI | 12 bp | `3M-5pgex-jan-2023` | No | Pending real-data validation |

## Barcode translation

For dual-oligo 3′ systems, the workflow performs translation in both directions:

1. GEX-to-feature translation for whitelist matching during guide extraction.
2. Feature-to-GEX translation after lane merge so guide and expression cell IDs
   use the same representation.

No translation is performed for 5′ single-oligo systems.

## Chemistry overrides

Individual fields can be overridden under `chemistry_overrides` without
defining a fully custom chemistry. Available keys are `af_chemistry`,
`whitelist`, `expected_ori`, `translation`, `translation_file`,
`geometry_override`, `ham_chemistry`, and `umi_len`.

## Custom chemistry

Use `custom` when the library does not match a supported preset:

```yaml
tenx_chemistry: custom
custom_chemistry:
  af_chemistry: "1{b[14]u[8]x:}2{r:}"
  whitelist: /path/to/whitelist.txt
  expected_ori: fw
  translation: false
  ham_chemistry: custom
  umi_len: 8

hash_matcher:
  custom_params:
    cb_start: 0
    cb_end: 14
    umi_start: 14
    umi_end: 22
    window_start: 8
    window_end: 28
    guide_len: 20
```

All required custom fields must be supplied explicitly.

## Unsupported 3′ v2 design

`3v2` is intentionally absent. It predates direct capture and uses indirect GBC
capture, so the guide sequence is not present in the direct-capture guide FASTQ
R2 layout expected by this workflow.
