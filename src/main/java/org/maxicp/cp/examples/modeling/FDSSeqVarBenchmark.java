/*
 * MaxiCP is under MIT License
 * Copyright (c)  2024 UCLouvain
 *
 */

package org.maxicp.cp.examples.modeling;

import org.maxicp.ModelDispatcher;
import org.maxicp.cp.examples.utils.TSPTWInstance;
import org.maxicp.cp.modeling.ConcreteCPModel;
import org.maxicp.modeling.SeqVar;
import org.maxicp.modeling.algebra.integer.IntExpression;
import org.maxicp.search.DFSearch;
import org.maxicp.search.FDSSeqVar;
import org.maxicp.search.SearchStatistics;
import org.maxicp.search.blackbox.BlackBoxSearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.OptionalInt;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.IntStream;

import static org.maxicp.modeling.Factory.*;
import static org.maxicp.search.Searches.*;

/**
 * Experiments of the FDSSeqVar paper (Tables 3 and 4).
 * <pre>
 * java FDSSeqVarBenchmark &lt;instance&gt; &lt;search&gt; &lt;timeLimitSec&gt; [repeat]
 * </pre>
 * {@code instance}: a Dumas TSPTW instance ({@code n20w20.001}) or a Cordeau DARP instance ({@code pr01}),
 * optionally reduced to its first {@code k} visits or requests ({@code n20w20.001-12}, {@code pr01-7}).
 * {@code search}: {@code nary} (N-ary minimum-detour DFS), {@code primal} ({@code Searches.primalSeqVar}),
 * {@code uniform} (FDSSeqVar with uniform slot density, no lazy pool), {@code pool} (uniform slot density +
 * lazy pool), {@code fds} ({@code Searches.fds}: slot-splitting estimator + lazy pool), {@code bb}
 * ({@code model.blackBoxSearch()}). A time limit of 0 runs until optimality is proven. The run is repeated
 * {@code repeat} times in the same JVM (to warm up the JIT) and one line is printed per run:
 * instance search limit nSol best tFirst tBest tEnd failures completed (times in ms).
 */
public class FDSSeqVarBenchmark {

    record Problem(ModelDispatcher model, SeqVar[] routes, IntExpression objective, int[][] dist, int[] nodes,
                   int scaling) {
    }

    public static void main(String[] args) {
        String instance = args[0];
        String search = args[1];
        int limitSec = Integer.parseInt(args[2]);
        int repeat = args.length > 3 ? Integer.parseInt(args[3]) : 1;
        for (int r = 0; r < repeat; r++) {
            Problem p = instance.startsWith("pr") ? darp(instance) : tsptw(instance);
            System.out.println(instance + " " + search + " " + limitSec + " " + run(p, search, limitSec));
        }
    }

    private static String run(Problem p, String search, int limitSec) {
        long limitMs = limitSec * 1000L;
        long start = System.currentTimeMillis();
        long[] sol = {0, -1, -1, -1}; // number of solutions, best objective, time of first and best solution
        SearchStatistics stats;
        if (search.equals("bb")) {
            BlackBoxSearch bb = p.model.blackBoxSearch();
            bb.onIncumbent((seq, obj) -> {
                long t = System.currentTimeMillis() - start;
                if (sol[0]++ == 0) sol[2] = t;
                sol[1] = obj;
                sol[3] = t;
            });
            stats = bb.optimize(p.model.minimize(p.objective), limitSec);
        } else {
            ConcreteCPModel cp = p.model.cpInstantiate();
            DFSearch dfs = cp.dfSearch(branching(p, search));
            dfs.onSolution(() -> {
                long t = System.currentTimeMillis() - start;
                if (sol[0]++ == 0) sol[2] = t;
                sol[1] = p.objective.min();
                sol[3] = t;
            });
            stats = dfs.optimize(cp.minimize(p.objective),
                    s -> limitMs > 0 && System.currentTimeMillis() - start >= limitMs);
        }
        long end = System.currentTimeMillis() - start;
        return String.format(Locale.ROOT, "%d %.2f %d %d %d %d %b", sol[0], (double) sol[1] / p.scaling, sol[2], sol[3], end,
                stats.numberOfFailures(), stats.isCompleted());
    }

    private static Supplier<Runnable[]> branching(Problem p, String search) {
        switch (search) {
            case "nary":
                // node with the fewest insertions, one branch per insertion by increasing detour
                Function<Integer, Runnable[]> branches = branchesInsertingNode(p.routes, p.dist).get();
                return () -> {
                    OptionalInt node = nodeSelector(p.routes, p.nodes, (seq, u) -> seq.nInsert(u));
                    return node.isEmpty() ? EMPTY : branches.apply(node.getAsInt());
                };
            case "primal":
                return primalSeqVar(p.routes, p.dist);
            case "uniform":
                return new FDSSeqVar(p.routes, p.dist).ablation(false, false);
            case "pool":
                return new FDSSeqVar(p.routes, p.dist).ablation(false, true);
            case "fds":
                return fds(p.routes, p.dist);
            default:
                throw new IllegalArgumentException("unknown search " + search);
        }
    }

    /**
     * TSPTW model of {@link org.maxicp.cp.examples.modeling.tsptw.TSPTWSeqVar}, on the depot and the first k
     * customers of the instance.
     */
    private static Problem tsptw(String name) {
        String[] parts = name.split("-");
        TSPTWInstance inst = new TSPTWInstance("data/TSPTW/Dumas/" + parts[0] + ".txt");
        int k = parts.length > 1 ? Integer.parseInt(parts[1]) : inst.n - 2;
        // depot, customers 1..k, end depot (duplicated by TSPTWInstance at index n - 1)
        int[] keep = IntStream.concat(IntStream.rangeClosed(0, k), IntStream.of(inst.n - 1)).toArray();
        int n = keep.length;
        int[][] dist = new int[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dist[i][j] = inst.distMatrix[keep[i]][keep[j]];
            }
        }
        ModelDispatcher model = makeModelDispatcher();
        SeqVar tour = model.seqVar(n, 0, n - 1);
        IntExpression[] time = new IntExpression[n];
        for (int i = 0; i < n; i++) {
            model.add(require(tour, i));
            time[i] = i == 0 ? model.intVar(0, 0) : model.intVar(inst.earliest[keep[i]], inst.latest[keep[i]]);
        }
        model.add(transitionTimes(tour, time, dist));
        IntExpression totalDist = model.intVar(0, org.maxicp.util.Arrays.max(dist) * n);
        model.add(distance(tour, dist, totalDist));
        return new Problem(model, new SeqVar[]{tour}, totalDist, dist, IntStream.range(1, n - 1).toArray(), 1);
    }

    /**
     * DARP model of {@link DARPSeqVar}, on the first k requests of the instance.
     */
    private static Problem darp(String name) {
        String[] parts = name.split("-");
        DARPSeqVar.Instance inst = new DARPSeqVar.Instance("data/DARP/Cordeau2003/" + parts[0], 100);
        if (parts.length > 1) {
            int k = Integer.parseInt(parts[1]);
            ArrayList<DARPSeqVar.NodeData> nodes = new ArrayList<>(inst.nodeData.subList(0, k)); // pickups
            nodes.addAll(inst.nodeData.subList(inst.nRequest, inst.nRequest + k)); // drops
            nodes.addAll(inst.nodeData.subList(2 * inst.nRequest, inst.nodeData.size())); // depots
            inst.nodeData = nodes;
            inst.nRequest = k;
            inst.distMatrix = new int[nodes.size()][nodes.size()];
            for (int i = 0; i < nodes.size(); i++) {
                for (int j = 0; j < nodes.size(); j++) {
                    inst.distMatrix[i][j] = nodes.get(i).distanceCeil(nodes.get(j));
                }
            }
        }
        int n = inst.nRequest * 2 + inst.nVehicle * 2;
        int[] start = IntStream.range(0, inst.nVehicle).map(i -> inst.nRequest * 2 + i).toArray();
        int[] end = Arrays.stream(start).map(i -> i + inst.nVehicle).toArray();
        ModelDispatcher model = makeModelDispatcher();
        SeqVar[] routes = new SeqVar[inst.nVehicle];
        IntExpression[] distance = new IntExpression[inst.nVehicle];
        for (int v = 0; v < inst.nVehicle; v++) {
            routes[v] = model.seqVar(n, start[v], end[v]);
            distance[v] = model.intVar(0, inst.get(end[v]).twEnd());
        }
        IntExpression[] time = new IntExpression[n];
        int[] duration = new int[n];
        for (int i = 0; i < n; i++) {
            time[i] = model.intVar(inst.get(i).twStart(), inst.get(i).twEnd());
            duration[i] = inst.get(i).duration();
        }
        int[] pickup = IntStream.range(0, inst.nRequest).toArray();
        int[] drop = IntStream.range(inst.nRequest, 2 * inst.nRequest).toArray();
        int[] load = IntStream.range(0, inst.nRequest).map(i -> inst.get(i).load()).toArray();
        for (int v = 0; v < inst.nVehicle; v++) {
            model.add(transitionTimes(routes[v], time, inst.distMatrix, duration));
            model.add(cumulative(routes[v], pickup, drop, load, inst.capacity));
            model.add(distance(routes[v], inst.distMatrix, distance[v]));
            model.add(le(time[end[v]], plus(time[start[v]], inst.maxRouteDuration)));
        }
        for (int r = 0; r < inst.nRequest; r++) {
            model.add(le(time[drop[r]], plus(time[pickup[r]], duration[pickup[r]] + inst.maxRideTime)));
            model.add(le(plus(time[pickup[r]], duration[pickup[r]] + inst.distMatrix[pickup[r]][drop[r]]),
                    time[drop[r]]));
        }
        for (int i = 0; i < n; i++) {
            int node = i;
            model.add(eq(sum(model.intVarArray(inst.nVehicle, v -> routes[v].isNodeRequired(node))), 1));
        }
        return new Problem(model, routes, sum(distance), inst.distMatrix,
                IntStream.range(0, 2 * inst.nRequest).toArray(), 100);
    }
}
