# gn-vep-backfill image

Image for the Ensembl 111 GRCh38 `vep.annotation` backfill. Built `FROM` the
gn-vep image that will read the cache, so annotations come from exactly that
VEP build.

## Build

Multi-arch is required, not optional. gn-vep is published for both
architectures, and the backfill runs across both: Pass A lands on the
`cbio-dev` pool, which is Graviton (`r6g`/`r8g` = arm64), while
`gn-vep-111-grch38` itself runs on `t3.large` (amd64).

```sh
docker buildx build --platform linux/amd64,linux/arm64 \
  -f Dockerfile.vep-backfill \
  -t genomenexus/gn-vep-backfill:111-grch38-1 --push .
```

Verify both platforms landed:

```sh
docker buildx imagetools inspect genomenexus/gn-vep-backfill:111-grch38-1
```

## What it adds

| | why |
|---|---|
| `python3` | the gn-vep image has **no python at all**; the converter needs it |
| `jq` | sets `_id` from VEP's `id` and flags `ref == alt` |
| `mongosh` | runs the covered `_id` scan against prod |
| `mongodb-database-tools` | `mongoimport` into the target collection |
| `/opt/backfill/*` | `hgvsg2ensembl.py`, `scan_vep_annotation_ids.js`, `export_vep_annotation_ids.sh` |

No pip step: the converter is stdlib-only by design, including its faidx
reader, so there is no `pysam`/`pyfaidx` to compile per architecture.

Inherited from the base image and deliberately not re-staged: the VEP 111 CLI
at `/opt/vep/src/ensembl-vep/vep`, and 13 GB of plugin data at `/plugin-data`
(PolyPhen_SIFT and AlphaMissense). Staging those on a volume would duplicate
~13 GB for no gain.

## Build-time checks

The build fails rather than deferring to a Job failure if `python3`, `jq`,
`mongosh` or `mongoimport` is missing or built for the wrong architecture, if
the VEP binary or either plugin database is absent, or if the converter's own
self-check does not pass:

```
convert('1:g.100A>C')      == '1 100 100 N/C + 1:g.100A>C'
minimise(100,101,'CT','TT') == (100,100,'C','T')
```

The second is the delins minimal-representation fix, which is what makes
`allele_string` match database mode.

## Runtime user

Reverts to `vep` (uid 999) as the base image does. EFS access points
substitute their own POSIX uid/gid for the client's, so writes to the cache and
work volumes succeed without `fsGroup` — verified in the cache loader, which
logged `write probe OK as 999:999` against a directory owned by `1000:1000`
with mode `700`.
