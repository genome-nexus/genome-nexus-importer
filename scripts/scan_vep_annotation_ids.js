// Stream every _id from vep.annotation as a covered index scan.
//
// A bare find({}, {_id:1}) plans as PROJECTION_SIMPLE <- COLLSCAN, which drags
// the whole collection through the WiredTiger cache. Verified on prod
// gn-mongo-grch38-ensembl95-blue: the ranged form below plans as
// LIMIT <- PROJECTION_COVERED <- IXSCAN(_id_), touching only the 4.51G _id
// index and never fetching a document. mongoexport cannot be used -- it has no
// --hint and table-scans by default on WiredTiger.
//
// env: BATCH (default 100000), SLEEP_MS (250), RESUME_AFTER (''), MAX_KEYS (0 = all)

(function () {
  const coll = db.getCollection("vep.annotation");
  const BATCH = parseInt(process.env.BATCH || "100000", 10);
  const SLEEP_MS = parseInt(process.env.SLEEP_MS || "250", 10);
  const RESUME_AFTER = process.env.RESUME_AFTER || "";
  const MAX_KEYS = parseInt(process.env.MAX_KEYS || "0", 10);

  // Fail loudly rather than silently table-scanning if the plan ever regresses.
  const probe = coll.find({ _id: { $gt: "" } }, { _id: 1 }).hint({ _id: 1 })
                    .sort({ _id: 1 }).limit(1).explain("queryPlanner");
  const plan = JSON.stringify(probe.queryPlanner.winningPlan);
  if (plan.indexOf("FETCH") !== -1 || plan.indexOf("COLLSCAN") !== -1) {
    throw new Error("refusing to scan: plan is not covered -> " + plan);
  }

  let last = RESUME_AFTER;
  let n = 0;
  while (true) {
    const q = (last === "") ? {} : { _id: { $gt: last } };
    const batch = coll.find(q, { _id: 1 }).hint({ _id: 1 })
                      .sort({ _id: 1 }).limit(BATCH).toArray();
    if (batch.length === 0) break;
    for (let i = 0; i < batch.length; i++) print(batch[i]._id);
    last = batch[batch.length - 1]._id;
    n += batch.length;
    if (MAX_KEYS > 0 && n >= MAX_KEYS) break;
    if (batch.length < BATCH) break;
    sleep(SLEEP_MS);
  }
})();
