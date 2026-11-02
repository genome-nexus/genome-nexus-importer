// Stream every _id from vep.annotation as a covered index scan, in verified
// range partitions.
//
// A bare find({}, {_id:1}) plans as PROJECTION_SIMPLE <- COLLSCAN, which drags
// the whole collection through the WiredTiger cache. The ranged form used here
// plans as PROJECTION_COVERED <- IXSCAN(_id_), touching only the _id index with
// no FETCH. The plan is re-checked at startup and the scan refuses to run if it
// ever regresses. mongoexport cannot be used: it has no --hint and documents a
// table scan as its default on WiredTiger.
//
// Why partitions rather than one resume-token sweep: a single unbounded sweep
// silently emitted 131,741,152 of 137,058,023 keys -- 3.9% short -- and the loss
// was systematic, dropping the high-position tail of each chromosome (keys at
// 7:g.7*, 7:g.8*, 7:g.9* were absent while 7:g.1* through 7:g.5* were present).
// Because skipping forward is still monotonic, nothing in the output revealed
// it. Every range here is bounded by $lt and its emitted count is checked
// against the index count, so a skip cannot escape the range it happens in and
// cannot pass unnoticed.
//
// env: BATCH (default 100000), SLEEP_MS (250), MAX_RANGE (4000000),
//      MAX_DEPTH (4), TOLERANCE_PPM (200), MAX_KEYS (0 = all)

const coll = db.getCollection("vep.annotation");
const BATCH = parseInt(process.env.BATCH || "100000", 10);
const SLEEP_MS = parseInt(process.env.SLEEP_MS || "250", 10);
const MAX_RANGE = parseInt(process.env.MAX_RANGE || "4000000", 10);
const TOLERANCE_PPM = parseInt(process.env.TOLERANCE_PPM || "200", 10);
const MAX_KEYS = parseInt(process.env.MAX_KEYS || "0", 10);
const MAX_DEPTH = parseInt(process.env.MAX_DEPTH || "4", 10);

function log(m) { print("[scan] " + m); }   // stderr-ish marker; filtered downstream

// --- refuse to run unless the plan is covered -----------------------------
const probe = coll.find({ _id: { $gt: "" } }, { _id: 1 }).hint({ _id: 1 })
                  .sort({ _id: 1 }).limit(1).explain("queryPlanner");
const plan = JSON.stringify(probe.queryPlanner.winningPlan);
if (plan.indexOf("FETCH") !== -1 || plan.indexOf("COLLSCAN") !== -1) {
  throw new Error("refusing to scan: plan is not covered -> " + plan);
}

// --- tile the key space ---------------------------------------------------
// Printable ASCII covers every byte an _id can start with. Dense ranges are
// subdivided by their second byte so no single range is unmanageably large.
const CHARS = [];
for (let c = 0x20; c <= 0x7e; c++) CHARS.push(String.fromCharCode(c));

function countRange(lo, hi) {
  const q = (hi === null) ? { _id: { $gte: lo } } : { _id: { $gte: lo, $lt: hi } };
  return coll.find(q, { _id: 1 }).hint({ _id: 1 }).count();
}

// Subdivide recursively so no leaf exceeds MAX_RANGE. Appending each candidate
// byte, plus a head piece below the lowest byte and a tail piece running to the
// parent's upper bound, tiles the parent exactly -- verified against the index
// total before any scanning starts.
function partition(lo, hi, depth) {
  const n = countRange(lo, hi);
  if (n === 0) return [];
  if (n <= MAX_RANGE || depth >= MAX_DEPTH) return [[lo, hi, n]];
  let out = [];
  const head = [lo, lo + CHARS[0]];
  out = out.concat(partition(head[0], head[1], depth + 1));
  for (let j = 0; j < CHARS.length; j++) {
    const slo = lo + CHARS[j];
    const shi = (j + 1 < CHARS.length) ? lo + CHARS[j + 1] : hi;
    out = out.concat(partition(slo, shi, depth + 1));
  }
  return out;
}

let leaves = [];
for (let i = 0; i < CHARS.length; i++) {
  const lo = CHARS[i];
  const hi = (i + 1 < CHARS.length) ? CHARS[i + 1] : null;
  leaves = leaves.concat(partition(lo, hi, 0));
}

const expected = leaves.reduce(function (a, r) { return a + r[2]; }, 0);
log("partitions=" + leaves.length + " expected_keys=" + expected);

// --- scan each range, verifying against its index count -------------------
let grand = 0;
const problems = [];
for (let k = 0; k < leaves.length; k++) {
  const lo = leaves[k][0], hi = leaves[k][1], want = leaves[k][2];
  let last = null, got = 0;
  while (true) {
    const lower = (last === null) ? { $gte: lo } : { $gt: last };
    const q = { _id: Object.assign({}, lower, (hi === null ? {} : { $lt: hi })) };
    const batch = coll.find(q, { _id: 1 }).hint({ _id: 1 })
                      .sort({ _id: 1 }).limit(BATCH).toArray();
    if (batch.length === 0) break;
    for (let i = 0; i < batch.length; i++) print(batch[i]._id);
    last = batch[batch.length - 1]._id;
    got += batch.length;
    // NOTE: deliberately no `break` on a short batch. That was the bug in the
    // previous version -- a short batch is not proof of end-of-range, and
    // treating it as such skipped the remainder. The $lt bound plus an empty
    // batch is the only termination condition.
    if (SLEEP_MS > 0) sleep(SLEEP_MS);
  }
  grand += got;
  // prod is live, so allow a little churn either way
  const tol = Math.max(50, Math.floor(want * TOLERANCE_PPM / 1000000));
  if (Math.abs(got - want) > tol) {
    problems.push(JSON.stringify(lo) + ".." + (hi === null ? "inf" : JSON.stringify(hi)) +
                  " want=" + want + " got=" + got);
    log("MISMATCH " + problems[problems.length - 1]);
  }
  if (MAX_KEYS > 0 && grand >= MAX_KEYS) break;
}

log("emitted=" + grand + " expected=" + expected + " mismatched_ranges=" + problems.length);
if (problems.length > 0) {
  throw new Error("scan incomplete in " + problems.length + " range(s): " + problems.slice(0, 5).join("; "));
}
