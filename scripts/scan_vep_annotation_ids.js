// Stream every _id from vep.annotation through a single covered cursor.
//
// One server-side aggregation cursor, one index pass, in index order. There are
// no range bounds, no resume tokens and no partition seams -- nothing to lose
// keys at. The plan is PROJECTION_COVERED <- IXSCAN with no FETCH and no
// COLLSCAN, so document pages are never touched; only the _id index is read.
//
// Two earlier designs both silently under-enumerated this collection:
//
//   a resume-token sweep           emitted 131,741,152 of 137,058,023  (-3.9%)
//   verified range partitions      emitted 130,913,776 of 137,058,023  (-4.5%)
//
// The sweep treated a short batch as end-of-range and skipped the remainder,
// dropping the high-position tail of each chromosome. The partitioned version
// replaced that with per-range count checks, which every range passed -- because
// the check compared each range against its own count and nothing verified that
// the ranges covered the key space, so a tiling gap was self-consistent and
// invisible. Both exited 0.
//
// mongoexport is not usable: it has no --hint and documents a table scan as its
// default on WiredTiger, which would drag all 273G of documents through a cache
// holding production's working set.
//
// env: BATCH_SIZE (default 10000), TOLERANCE_PPM (200), MAX_KEYS (0 = all)

const coll = db.getCollection("vep.annotation");
const BATCH_SIZE = parseInt(process.env.BATCH_SIZE || "10000", 10);
const TOLERANCE_PPM = parseInt(process.env.TOLERANCE_PPM || "200", 10);
const MAX_KEYS = parseInt(process.env.MAX_KEYS || "0", 10);

// Diagnostics are prefixed so the downstream awk filter can forward them to
// stderr rather than mistaking them for keys.
function log(m) { print("[scan] " + m); }

const PIPELINE = [{ $project: { _id: 1 } }];
const OPTS = { hint: { _id: 1 }, batchSize: BATCH_SIZE, allowDiskUse: false };

// --- refuse to run unless the plan is covered -----------------------------
const plan = JSON.stringify(coll.aggregate(PIPELINE, OPTS).explain("queryPlanner"));
if (plan.indexOf("FETCH") !== -1 || plan.indexOf("COLLSCAN") !== -1) {
  throw new Error("refusing to scan: plan is not covered -> " + plan);
}

// Authoritative total, taken before the pass. Prod is live, so a small drift
// either way is expected; a real shortfall is orders of magnitude larger.
const expected = coll.find({ _id: { $gte: "" } }, { _id: 1 }).hint({ _id: 1 }).count();
log("expected_keys=" + expected + " batch_size=" + BATCH_SIZE);

// --- the single pass ------------------------------------------------------
const cursor = coll.aggregate(PIPELINE, OPTS);
let n = 0;
while (cursor.hasNext()) {
  print(cursor.next()._id);
  n++;
  if (MAX_KEYS > 0 && n >= MAX_KEYS) break;
}

log("emitted=" + n + " expected=" + expected);

// --- verify the whole pass, not pieces of it ------------------------------
if (MAX_KEYS === 0) {
  const tol = Math.max(1000, Math.floor(expected * TOLERANCE_PPM / 1000000));
  if (expected - n > tol) {
    throw new Error("scan incomplete: emitted " + n + " of " + expected +
                    " (short by " + (expected - n) + ", tolerance " + tol + ")");
  }
}
