# gn-vep-backfill image

Image for the Ensembl 111 GRCh38 `vep.annotation` backfill. Built `FROM` the
gn-vep image that will read the cache, so annotations come from exactly that
VEP build.

## Build

Published as `genomenexus/gn-vep-backfill:111-grch38-1`, multi-arch
(`linux/amd64` `sha256:70f9ff62`, `linux/arm64` `sha256:7325db77`).

```sh
docker buildx build --platform linux/amd64,linux/arm64 \
  -f Dockerfile.vep-backfill \
  -t genomenexus/gn-vep-backfill:111-grch38-1 --push .
```

Verify:

```sh
docker buildx imagetools inspect genomenexus/gn-vep-backfill:111-grch38-1
```

The Docker Hub web API lags behind a push by a few minutes; the registry is
authoritative if you need to confirm sooner:

```sh
TOK=$(curl -s "https://auth.docker.io/token?service=registry.docker.io&scope=repository:genomenexus/gn-vep-backfill:pull" | jq -r .token)
curl -s -H "Authorization: Bearer $TOK" \
  https://registry-1.docker.io/v2/genomenexus/gn-vep-backfill/tags/list
```

### If you rebuild

Two things make it much faster.

**Tag inside the `gn-vep` repository** rather than a separate one. A registry
only deduplicates blobs within the *same* repository name, so pushing to a new
repo re-uploads the entire 12.4 GB base per platform — measured at over half an
hour. Tagging as `gn-vep:111-grch38-backfill-N` makes every base blob return 200
from the registry's `HEAD` check, leaving only a few hundred MB to transfer.

**`--platform linux/arm64` alone is sufficient.** The backfill Jobs target the
`cbio-dev` pool, which is Graviton (`r6g`/`r8g`). amd64 is built only so the
image remains usable on the amd64 pools; nothing in the backfill needs it.

Build locally rather than in CI. `.github/workflows/docker-image-push.yml` notes
that GitHub Actions cannot build and push images this size, and its `type=gha`
cache would exceed the 10 GB limit.

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
