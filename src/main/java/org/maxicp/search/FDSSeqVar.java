/*
 * MaxiCP is under MIT License
 * Copyright (c)  2024 UCLouvain
 */

package org.maxicp.search;

import org.maxicp.cp.engine.core.CPSeqVar;
import org.maxicp.cp.engine.core.CPSolver;
import org.maxicp.cp.modeling.ConcreteCPModel;
import org.maxicp.modeling.Factory;
import org.maxicp.modeling.ModelProxy;
import org.maxicp.modeling.SeqVar;
import org.maxicp.modeling.concrete.ConcreteModel;
import org.maxicp.state.StateManager;
import org.maxicp.util.TriFunction;
import org.maxicp.util.exception.InconsistencyException;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.function.Supplier;

import static org.maxicp.modeling.algebra.sequence.SeqStatus.*;
import static org.maxicp.search.Searches.EMPTY;
import static org.maxicp.search.Searches.branch;

/**
 * Failure-Directed Search (FDS) adapted for Sequence Variables ({@link SeqVar} and {@link CPSeqVar}).
 * <p>
 * This search strategy adapts the principles of Failure-Directed Search (Vilím, Laborie, Shaw, CPAIOR 2015)
 * to sequence variables using binary insertion decisions:
 * <ul>
 *   <li><b>Positive branch (+):</b> {@code insert(seq, pred, node)} — inserts {@code node} inside the current
 *       member interval {@code (pred, succ)}, where {@code succ = seq.memberAfter(pred)}.</li>
 *   <li><b>Negative branch (-):</b> {@code notBetween(seq, pred, node, succ)} — forbids {@code node} from ever
 *       being inserted inside the interval {@code (pred, succ)}.</li>
 * </ul>
 * <p>
 * Key features of Sequence FDS:
 * <ol>
 *   <li><b>Exact Binary Partition &amp; Active Refutation Propagation:</b> Unlike N-ary branching over all feasible
 *       insertions of a node, refuting the positive branch explicitly posts {@code notBetween(seq, pred, node, succ)}
 *       and runs {@code fixPoint()} on the negative branch. This triggers {@code propagateOnInsertRemoved} in sequence
 *       constraints ({@code Cumulative}, {@code Precedence}), automatically excludes {@code node} from {@code seq}
 *       when its insertion count drops to 0, and automatically inserts {@code node} when it is required and only
 *       1 valid insertion remains.</li>
 *   <li><b>Dynamic Choice Activation &amp; Transition Sharing:</b> Choices {@code (pred, node)} become active
 *       dynamically as {@code pred} becomes a member of a sequence ({@code seq.hasInsert(pred, node)}). Indexing
 *       ratings by transition {@code (pred, node)} shares learned conflict ratings across symmetrical vehicles
 *       while keeping vehicle-specific depot transitions distinct.</li>
 *   <li><b>Combinatorial Space Contraction via Sequential Slot-Splitting ({@code R}):</b> The exact number of valid
 *       sequence completions of a sequence domain with {@code K} open member intervals (slots) is
 *       {@code sum_f prod_{j=1}^K (r_j(f))!}, which is #P-complete to compute in general. {@code FDSSeqVar}
 *       approximates the number of options {@code choices(u)} and the slot density {@code rho_seq(u)} of each
 *       uninserted node in {@code O(|U| * K)} time with the Sequential Slot-Splitting Estimator, which is exact on
 *       laminar domains without possible nodes. The slots of all sequences form a single pool, and a node counts as
 *       required if it is required in one of the sequences where it is insertable.</li>
 *   <li><b>Predecessor Conflict Blame:</b> When {@code insert(seq, pred, node)} fails, conflict blame is also
 *       propagated to {@code predRating[pred]}, accelerating failure detection on coupled pairs (such as
 *       pickup-and-delivery requests in DARP).</li>
 *   <li><b>Short-Detour Lazy Choice Pool:</b> When a distance matrix or detour cost function
 *       {@code (pred, node, succ) -> cost} is supplied, {@code FDSSeqVar} initially only admits choice points
 *       {@code (pred, node)} with a short detour ({@code <= tau_pool * avgEdgeCost}), the minimum-detour insertion of
 *       each node on each sequence, and the transitions already selected for branching ({@code inPool[pred][node]}).</li>
 * </ol>
 *
 * @see FDS
 * @see CPSeqVar
 * @see SeqVar
 */
public class FDSSeqVar implements Supplier<Runnable[]> {

    private static final double DEFAULT_ALPHA = 0.95;
    private static final double DEFAULT_DETOUR_WEIGHT = 0.20;
    private static final double DEFAULT_INITIAL_POOL_TAU = 0.75;

    private final SeqVar[] seqVars;
    private final ModelProxy model;
    private final double alpha;
    private double detourWeight;
    private final double initialPoolTau;
    /** Ablation switches ({@link #ablation}): slot-splitting weight updates and the short-detour lazy pool. */
    private boolean slotSplitting = true;
    private boolean lazyPool = true;
    private final boolean autoDiscoverDistMatrix;
    private int[][] distMatrix;
    private final TriFunction<Integer, Integer, Integer, Integer> detourCost;
    private final int nTotalNodes;

    /** Positive branch ratings indexed by [pred][node]. */
    private final double[][] ratingPos;
    /** Negative branch ratings indexed by [pred][node]. */
    private final double[][] ratingNeg;
    /** Predecessor conflict rating indexed by [pred], tracking downstream failures after pred is inserted. */
    private final double[] predRating;
    /** Lazy choice pool indicator indexed by [pred][node] when detour guidance is active. */
    private final boolean[][] inPool;
    /** Running average of local ratings per search depth: [0] = mean, [1] = count. */
    private final Map<Integer, double[]> avgRatingByDepth;

    /** Zero-allocation reusable buffers for sequential slot-splitting snapshots and candidate pool selection. */
    private final int[] insertions;
    private final int[] activeSlots;
    private final double[] rhoBefore;
    private final double[] choicesBefore;
    private final double[] rhoAfter;
    private final double[] choicesAfter;
    private final double[] slotWeights;
    private final int[] nodeInsertCounts;
    private final int[] bucketHeads;
    private final int[] bucketNext;
    private final int[] candSeq;
    private final int[] candPred;
    private final int[] candDetour;
    private final boolean[] candAdmitted;
    private double avgEdgeCost = 1.0;

    /** Cached concrete model and raw CP variables for direct zero-overhead access. */
    private ConcreteModel cachedConcreteModel;
    private CPSolver cpSolver;
    private CPSeqVar[] cpSeqVars;
    private StateManager sm;
    private int baseLevel;

    /**
     * Creates a generic Failure-Directed Search for the given sequence variables.
     *
     * @param seqVars the sequence variables to decide
     */
    public FDSSeqVar(SeqVar... seqVars) {
        this(DEFAULT_ALPHA, 0.0, seqVars, null, null, true);
    }

    /**
     * Creates a generic Failure-Directed Search for the given sequence variables with custom decay factor.
     *
     * @param alpha   exponential decay factor for rating updates (typical: 0.9 to 0.99)
     * @param seqVars the sequence variables to decide
     */
    public FDSSeqVar(double alpha, SeqVar... seqVars) {
        this(alpha, 0.0, seqVars, null, null, true);
    }

    /**
     * Creates a Failure-Directed Search for the given sequence variables guided by a distance matrix.
     * The detour cost of inserting {@code node} between {@code pred} and {@code succ} is computed as
     * {@code distMatrix[pred][node] + distMatrix[node][succ] - distMatrix[pred][succ]}.
     *
     * @param seqVars    the sequence variables to decide
     * @param distMatrix the distance matrix between nodes
     */
    public FDSSeqVar(SeqVar[] seqVars, int[][] distMatrix) {
        this(DEFAULT_ALPHA, DEFAULT_DETOUR_WEIGHT, seqVars, distMatrix, null, false);
    }

    /**
     * Creates a Failure-Directed Search for the given sequence variables guided by a custom detour cost function.
     *
     * @param seqVars    the sequence variables to decide
     * @param detourCost function computing the detour cost for {@code (pred, node, succ)}
     */
    public FDSSeqVar(SeqVar[] seqVars, TriFunction<Integer, Integer, Integer, Integer> detourCost) {
        this(DEFAULT_ALPHA, DEFAULT_DETOUR_WEIGHT, seqVars, null, detourCost, false);
    }

    /**
     * Creates a Failure-Directed Search for the given sequence variables with custom parameters and distance matrix.
     *
     * @param alpha        exponential decay factor for FDS rating updates (typical: 0.9 to 0.99)
     * @param detourWeight weight of the normalized detour cost relative to the FDS rating (typical: 0.05 to 0.25)
     * @param seqVars      the sequence variables to decide
     * @param distMatrix   the distance matrix between nodes (may be null)
     */
    public FDSSeqVar(double alpha, double detourWeight, SeqVar[] seqVars, int[][] distMatrix) {
        this(alpha, detourWeight, seqVars, distMatrix, null, false);
    }

    /**
     * Creates a Failure-Directed Search for the given sequence variables with custom parameters and detour function.
     *
     * @param alpha        exponential decay factor for FDS rating updates (typical: 0.9 to 0.99)
     * @param detourWeight weight of the normalized detour cost relative to the FDS rating (typical: 0.05 to 0.25)
     * @param seqVars      the sequence variables to decide
     * @param detourCost   function computing the detour cost for {@code (pred, node, succ)} (may be null)
     */
    public FDSSeqVar(double alpha, double detourWeight, SeqVar[] seqVars,
                     TriFunction<Integer, Integer, Integer, Integer> detourCost) {
        this(alpha, detourWeight, seqVars, null, detourCost, false);
    }

    private FDSSeqVar(double alpha, double detourWeight, SeqVar[] seqVars,
                      int[][] distMatrix, TriFunction<Integer, Integer, Integer, Integer> detourCost,
                      boolean autoDiscoverDistMatrix) {
        if (seqVars == null || seqVars.length == 0) {
            throw new IllegalArgumentException("FDSSeqVar requires at least one sequence variable");
        }
        for (SeqVar sv : seqVars) {
            Objects.requireNonNull(sv, "Sequence variable cannot be null");
        }
        if (alpha <= 0.0 || alpha >= 1.0) {
            throw new IllegalArgumentException("alpha must be in (0, 1)");
        }
        this.seqVars = seqVars;
        this.alpha = alpha;
        this.detourWeight = Math.max(0.0, detourWeight);
        this.initialPoolTau = DEFAULT_INITIAL_POOL_TAU;
        this.autoDiscoverDistMatrix = autoDiscoverDistMatrix;
        this.distMatrix = distMatrix;
        this.detourCost = detourCost;
        this.model = resolveModelProxy(seqVars);

        int maxNodes = 0;
        for (SeqVar sv : seqVars) {
            maxNodes = Math.max(maxNodes, sv.nNode());
        }
        this.nTotalNodes = maxNodes;

        this.ratingPos = new double[nTotalNodes][nTotalNodes];
        this.ratingNeg = new double[nTotalNodes][nTotalNodes];
        for (int i = 0; i < nTotalNodes; i++) {
            Arrays.fill(ratingPos[i], 1.0);
            Arrays.fill(ratingNeg[i], 1.0);
        }
        this.predRating = new double[nTotalNodes];
        Arrays.fill(predRating, 1.0);
        this.inPool = new boolean[nTotalNodes][nTotalNodes];
        this.avgRatingByDepth = new HashMap<>();

        int maxPossibleIns = seqVars.length * nTotalNodes + 2;
        this.insertions = new int[nTotalNodes];
        this.activeSlots = new int[maxPossibleIns];
        this.rhoBefore = new double[nTotalNodes];
        this.choicesBefore = new double[nTotalNodes];
        this.rhoAfter = new double[nTotalNodes];
        this.choicesAfter = new double[nTotalNodes];
        this.slotWeights = new double[seqVars.length * nTotalNodes];
        this.nodeInsertCounts = new int[nTotalNodes];
        this.bucketHeads = new int[maxPossibleIns];
        this.bucketNext = new int[nTotalNodes];
        this.candSeq = new int[maxPossibleIns];
        this.candPred = new int[maxPossibleIns];
        this.candDetour = new int[maxPossibleIns];
        this.candAdmitted = new boolean[maxPossibleIns];

        if (distMatrix != null) {
            updateAvgEdgeCostFromMatrix(distMatrix);
        } else if (detourCost != null) {
            long sum = 0;
            int cnt = 0;
            int limit = Math.min(nTotalNodes, 32);
            for (int i = 0; i < limit; i++) {
                for (int j = 0; j < limit; j++) {
                    if (i == j) continue;
                    int c = Math.abs(detourCost.apply(i, j, i));
                    if (c > 0) {
                        sum += c;
                        cnt++;
                    }
                }
            }
            if (cnt > 0) {
                // detourCost(i, j, i) = d[i][j] + d[j][i]: halve it to get the mean transition distance
                this.avgEdgeCost = (double) sum / (2.0 * cnt);
            }
        }

        if (this.model == null || this.model.isConcrete()) {
            ensureConcreteCache();
        }
    }

    private void updateAvgEdgeCostFromMatrix(int[][] matrix) {
        long sum = 0;
        int cnt = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (i != j && matrix[i][j] > 0) {
                    sum += matrix[i][j];
                    cnt++;
                }
            }
        }
        if (cnt > 0) {
            this.avgEdgeCost = (double) sum / cnt;
        }
    }

    private static ModelProxy resolveModelProxy(SeqVar[] seqVars) {
        boolean allRaw = true;
        for (SeqVar sv : seqVars) {
            if (!(sv instanceof CPSeqVar)) {
                allRaw = false;
                break;
            }
        }
        if (allRaw) {
            CPSolver solver = ((CPSeqVar) seqVars[0]).getSolver();
            for (int i = 1; i < seqVars.length; i++) {
                if (((CPSeqVar) seqVars[i]).getSolver() != solver) {
                    throw new IllegalArgumentException("All sequence variables must belong to the same solver");
                }
            }
            return null;
        }
        ModelProxy proxy = null;
        for (SeqVar sv : seqVars) {
            if (proxy == null) {
                proxy = sv.getModelProxy();
            } else if (sv.getModelProxy() != null && proxy != sv.getModelProxy()) {
                throw new IllegalArgumentException("All sequence variables must belong to the same model");
            }
        }
        return proxy;
    }

    private void ensureConcreteCache() {
        if (model == null) {
            if (cpSolver == null) {
                cpSeqVars = new CPSeqVar[seqVars.length];
                for (int i = 0; i < seqVars.length; i++) {
                    cpSeqVars[i] = (CPSeqVar) seqVars[i];
                }
                cpSolver = cpSeqVars[0].getSolver();
                sm = cpSolver.getStateManager();
                baseLevel = sm.getLevel();
                tryDiscoverDistMatrix();
            }
            return;
        }
        ConcreteModel cm = model.getConcreteModel();
        if (cm != cachedConcreteModel) {
            cachedConcreteModel = cm;
            sm = cm.getStateManager();
            baseLevel = sm.getLevel();
            if (cm instanceof ConcreteCPModel cpm) {
                cpSolver = cpm.solver;
                cpSeqVars = new CPSeqVar[seqVars.length];
                for (int i = 0; i < seqVars.length; i++) {
                    cpSeqVars[i] = cpm.getCPVar(seqVars[i]);
                }
                tryDiscoverDistMatrix();
            } else {
                cpSolver = null;
                cpSeqVars = null;
            }
        }
    }

    private void tryDiscoverDistMatrix() {
        if (!autoDiscoverDistMatrix || distMatrix != null || detourCost != null || cpSeqVars == null) {
            return;
        }
        for (CPSeqVar cpSeq : cpSeqVars) {
            if (cpSeq != null && cpSeq.distanceMatrix() != null) {
                this.distMatrix = cpSeq.distanceMatrix();
                this.detourWeight = DEFAULT_DETOUR_WEIGHT;
                updateAvgEdgeCostFromMatrix(this.distMatrix);
                return;
            }
        }
    }

    /**
     * Resets the depth base level to the current solver level while keeping learned choice ratings.
     */
    public void resetLevel() {
        ensureConcreteCache();
        if (sm != null) {
            this.baseLevel = sm.getLevel();
        }
    }

    /**
     * Disables components for ablation experiments. With {@code slotSplitting = false}, the reduction ratio
     * {@code R} uses the uniform slot density {@code nInsert(u) / K} (the slot weights are never updated).
     * With {@code lazyPool = false}, every feasible insertion is scored instead of the short-detour lazy pool.
     *
     * @return this search, for chaining
     */
    public FDSSeqVar ablation(boolean slotSplitting, boolean lazyPool) {
        this.slotSplitting = slotSplitting;
        this.lazyPool = lazyPool;
        return this;
    }

    private int currentDepth() {
        return sm == null ? 0 : Math.max(0, sm.getLevel() - baseLevel);
    }

    private double[] getAvgRating(int depth) {
        return avgRatingByDepth.computeIfAbsent(depth, k -> new double[]{1.0, 0.0});
    }

    private SeqVar seqAt(int idx) {
        return cpSeqVars != null ? cpSeqVars[idx] : seqVars[idx];
    }

    /**
     * Computes the number of options {@code choices(u)} and the sequential slot density
     * {@code rho_seq(u) = W_u / W} of each uninserted node {@code u} with the Sequential Slot-Splitting
     * Estimator in {@code O(|U| * K)} time.
     * <p>
     * Nodes are processed by non-decreasing {@code nInsert(u)} (ties by node index). Node {@code u} sees
     * {@code W_u = sum_{j in J_u} w_j} expected gaps: {@code choices(u) = W_u} if {@code u} is required and
     * {@code 1 + W_u} if it is possible (it may also be excluded). It is inserted with probability
     * {@code p_u = 1} (required) or {@code W_u / (1 + W_u)} (possible), so {@code w_j <- w_j * (1 + p_u / W_u)}
     * for {@code j in J_u} and {@code W <- W + p_u}.
     */
    private void snapshotSequentialSlotSplitting(double[] destRho, double[] destChoices) {
        int nSeq = seqVars.length;
        int totalSlots = 0;
        for (int s = 0; s < nSeq; s++) {
            totalSlots += Math.max(1, seqAt(s).nNode(MEMBER) - 1);
        }
        totalSlots = Math.max(1, totalSlots);

        int maxIns = 0;
        for (int u = 0; u < nTotalNodes; u++) {
            int ins = 0;
            for (int s = 0; s < nSeq; s++) {
                SeqVar r = seqAt(s);
                if (u < r.nNode() && r.isNode(u, INSERTABLE)) {
                    ins += r.nInsert(u);
                }
            }
            nodeInsertCounts[u] = ins;
            destRho[u] = 0.0;
            destChoices[u] = 0.0;
            if (ins > maxIns) {
                maxIns = ins;
            }
        }
        if (maxIns == 0) {
            return;
        }

        Arrays.fill(bucketHeads, 0, maxIns + 1, -1);
        for (int u = nTotalNodes - 1; u >= 0; u--) {
            int ins = nodeInsertCounts[u];
            if (ins > 0) {
                bucketNext[u] = bucketHeads[ins];
                bucketHeads[ins] = u;
            }
        }

        Arrays.fill(slotWeights, 1.0);
        double totalWeight = totalSlots;

        for (int d = 1; d <= maxIns; d++) {
            for (int u = bucketHeads[d]; u != -1; u = bucketNext[u]) {
                double wU = 0.0;
                boolean required = false;
                int nSlotsU = 0;
                for (int s = 0; s < nSeq; s++) {
                    SeqVar r = seqAt(s);
                    if (u >= r.nNode() || !r.isNode(u, INSERTABLE)) continue;
                    required |= r.isNode(u, REQUIRED);
                    int nIns = r.fillInsert(u, insertions);
                    int base = s * nTotalNodes;
                    for (int j = 0; j < nIns; j++) {
                        int flatSlot = base + insertions[j];
                        activeSlots[nSlotsU++] = flatSlot;
                        wU += slotWeights[flatSlot];
                    }
                }
                if (wU <= 0.0) continue;
                destChoices[u] = required ? wU : 1.0 + wU;
                destRho[u] = wU / totalWeight;
                if (!slotSplitting) continue; // uniform slot density: all weights stay 1
                double p = required ? 1.0 : wU / (1.0 + wU);
                double factor = 1.0 + p / wU;
                for (int k = 0; k < nSlotsU; k++) {
                    slotWeights[activeSlots[k]] *= factor;
                }
                totalWeight += p;
            }
        }
    }

    /**
     * Computes the search space reduction ratio {@code R} in {@code [0, 1]} after a branch and fixpoint
     * using the Sequential Slot-Splitting Estimator.
     */
    private double computeSequentialSlotReduction() {
        snapshotSequentialSlotSplitting(rhoAfter, choicesAfter);
        double ratio = 1.0;
        for (int u = 0; u < nTotalNodes; u++) {
            double cBefore = choicesBefore[u];
            if (cBefore > 1.0 + 1e-9) {
                double cAfter = choicesAfter[u];
                if (cAfter <= 0.0) {
                    ratio *= 1.0 / cBefore;
                } else {
                    double nodeRatio = rhoAfter[u] / rhoBefore[u];
                    if (nodeRatio < 1.0) {
                        ratio *= nodeRatio;
                    }
                }
            }
        }
        return Math.min(ratio, 1.0);
    }

    private void updateRating(int pred, int node, boolean isPositive, boolean failed, double R) {
        double localRating = failed ? 0.0 : 1.0 + R;
        int d = currentDepth();
        double[] avg = getAvgRating(d);
        avg[1] += 1.0;
        avg[0] += (localRating - avg[0]) / avg[1];
        double avgVal = avg[0] > 1e-6 ? avg[0] : 1.0;
        double norm = localRating / avgVal;

        if (isPositive) {
            ratingPos[pred][node] = alpha * ratingPos[pred][node] + (1.0 - alpha) * norm;
        } else {
            ratingNeg[pred][node] = alpha * ratingNeg[pred][node] + (1.0 - alpha) * norm;
        }
        predRating[pred] = alpha * predRating[pred] + (1.0 - alpha) * norm;
    }

    private int evaluateDetour(SeqVar route, int pred, int node, int succ) {
        if (distMatrix != null) {
            return distMatrix[pred][node] + distMatrix[node][succ] - distMatrix[pred][succ];
        }
        if (detourCost != null) {
            return detourCost.apply(pred, node, succ);
        }
        return route.nPred(node) + route.nSucc(node);
    }

    private boolean hasDetourGuidance() {
        return (distMatrix != null || detourCost != null) && detourWeight > 0.0;
    }

    @Override
    public Runnable[] get() {
        ensureConcreteCache();

        int nSeq = seqVars.length;
        boolean anyUnfixed = false;
        boolean anyRequiredUninserted = false;

        for (int u = 0; u < nTotalNodes; u++) {
            for (int s = 0; s < nSeq; s++) {
                SeqVar r = seqAt(s);
                if (u >= r.nNode()) continue;
                if (r.isNode(u, REQUIRED) && !r.isNode(u, MEMBER) && r.nInsert(u) == 0) {
                    return new Runnable[]{() -> {
                        throw InconsistencyException.INCONSISTENCY;
                    }};
                }
                if (r.isNode(u, INSERTABLE)) {
                    int nIns = r.nInsert(u);
                    if (nIns > 0) {
                        anyUnfixed = true;
                        if (r.isNode(u, REQUIRED)) {
                            anyRequiredUninserted = true;
                        }
                    }
                }
            }
        }
        if (!anyUnfixed) {
            return EMPTY;
        }

        int bestSeqIdx = -1;
        int bestPred = -1;
        int bestNode = -1;
        double bestScore = Double.MAX_VALUE;
        int bestTie = Integer.MAX_VALUE;
        boolean useDetour = hasDetourGuidance();
        boolean usePool = useDetour && lazyPool;
        double fdsWeight = 1.0 - detourWeight;
        double shortThreshold = initialPoolTau * avgEdgeCost;

        for (int u = 0; u < nTotalNodes; u++) {
            int totalIns = 0;
            boolean isReq = false;
            boolean isIns = false;
            for (int s = 0; s < nSeq; s++) {
                SeqVar r = seqAt(s);
                if (u >= r.nNode()) continue;
                if (r.isNode(u, INSERTABLE)) {
                    int nIns = r.nInsert(u);
                    if (nIns > 0) {
                        isIns = true;
                        totalIns += nIns;
                        if (r.isNode(u, REQUIRED)) {
                            isReq = true;
                        }
                    }
                }
            }
            if (!isIns) continue;
            if (anyRequiredUninserted && !isReq) continue;

            int nCand = 0;
            for (int s = 0; s < nSeq; s++) {
                SeqVar r = seqAt(s);
                if (u >= r.nNode() || !r.isNode(u, INSERTABLE)) continue;
                if (isReq && !r.isNode(u, REQUIRED)) continue;
                int nIns = r.fillInsert(u, insertions);
                int argMin = nCand; // minimum-detour insertion of u on this sequence, always admitted
                for (int j = 0; j < nIns; j++) {
                    int pred = insertions[j];
                    int succ = r.memberAfter(pred);
                    int detour = evaluateDetour(r, pred, u, succ);
                    candSeq[nCand] = s;
                    candPred[nCand] = pred;
                    candDetour[nCand] = detour;
                    candAdmitted[nCand] = !usePool || inPool[pred][u] || detour <= shortThreshold;
                    if (detour < candDetour[argMin]) {
                        argMin = nCand;
                    }
                    nCand++;
                }
                if (nIns > 0) {
                    candAdmitted[argMin] = true;
                }
            }

            for (int c = 0; c < nCand; c++) {
                if (!candAdmitted[c]) continue;
                int pred = candPred[c];

                int s = candSeq[c];
                int tie = candDetour[c];
                double rawFds = (ratingPos[pred][u] + ratingNeg[pred][u]) * (0.9 + 0.1 * predRating[u]);
                double cRating;
                if (useDetour) {
                    double detourNorm = tie / avgEdgeCost;
                    cRating = totalIns * (fdsWeight * rawFds + detourWeight * detourNorm);
                } else {
                    cRating = totalIns * rawFds;
                }
                if (cRating < bestScore - 1e-9 || (Math.abs(cRating - bestScore) <= 1e-9 && tie < bestTie)) {
                    bestScore = cRating;
                    bestTie = tie;
                    bestSeqIdx = s;
                    bestPred = pred;
                    bestNode = u;
                }
            }
        }

        if (bestSeqIdx == -1) {
            return new Runnable[]{() -> {
                throw InconsistencyException.INCONSISTENCY;
            }};
        }

        final int seqIdx = bestSeqIdx;
        final int pred = bestPred;
        final int node = bestNode;
        final int succ = seqAt(seqIdx).memberAfter(pred);
        inPool[pred][node] = true; // a transition selected for branching stays in the pool
        boolean posFirst = useDetour || (ratingPos[pred][node] <= ratingNeg[pred][node]);
        Runnable posBranch = makeBranch(seqIdx, pred, node, succ, true);
        Runnable negBranch = makeBranch(seqIdx, pred, node, succ, false);
        return posFirst ? branch(posBranch, negBranch) : branch(negBranch, posBranch);
    }

    private Runnable makeBranch(int seqIdx, int pred, int node, int succ, boolean isPositive) {
        return () -> {
            snapshotSequentialSlotSplitting(rhoBefore, choicesBefore);
            try {
                if (cpSeqVars != null) {
                    CPSeqVar route = cpSeqVars[seqIdx];
                    if (isPositive) {
                        if (!route.hasInsert(pred, node)) {
                            throw InconsistencyException.INCONSISTENCY;
                        }
                        route.insert(pred, node);
                    } else {
                        if (!route.isNode(pred, MEMBER) || !route.isNode(succ, MEMBER)) {
                            throw InconsistencyException.INCONSISTENCY;
                        }
                        route.notBetween(pred, node, succ);
                    }
                    cpSolver.fixPoint();
                } else {
                    SeqVar route = seqVars[seqIdx];
                    if (isPositive) {
                        if (!route.hasInsert(pred, node)) {
                            throw InconsistencyException.INCONSISTENCY;
                        }
                        model.add(Factory.insert(route, pred, node));
                    } else {
                        if (!route.isNode(pred, MEMBER) || !route.isNode(succ, MEMBER)) {
                            throw InconsistencyException.INCONSISTENCY;
                        }
                        model.add(Factory.notBetween(route, pred, node, succ));
                    }
                }
            } catch (InconsistencyException e) {
                updateRating(pred, node, isPositive, true, 0.0);
                throw e;
            }
            double R = computeSequentialSlotReduction();
            updateRating(pred, node, isPositive, false, R);
        };
    }
}
