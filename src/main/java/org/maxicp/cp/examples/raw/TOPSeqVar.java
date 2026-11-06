package org.maxicp.cp.examples.raw;

import org.maxicp.cp.engine.constraints.Sum;
import org.maxicp.cp.engine.constraints.seqvar.Distance;
import org.maxicp.cp.engine.core.CPIntVar;
import org.maxicp.cp.engine.core.CPSeqVar;
import org.maxicp.cp.engine.core.CPSolver;
import org.maxicp.cp.examples.utils.TOPInstance;
import org.maxicp.cp.examples.utils.TopFormat;
import org.maxicp.search.DFSearch;
import org.maxicp.search.Objective;
import org.maxicp.search.SearchStatistics;

import java.io.File;
import java.util.*;

import static org.maxicp.cp.CPFactory.*;
import static org.maxicp.modeling.algebra.sequence.SeqStatus.*;
import static org.maxicp.search.Searches.*;

public class TOPSeqVar {

    // number of destroyed customers, in % of the customers (Chaigneau et al. 2025, Table 2)
    public static final int MIN_CUSTOMER_DESTROYED = 1;
    public static final int MAX_CUSTOMER_DESTROYED = 20;

    public static final int MIN_SEQ_CUSTOMER_DESTROYED = 1;
    public static final int MAX_SEQ_CUSTOMER_DESTROYED = 20;

    public static final int LNS_RUNNING_TIME = 95*1000; // in ms
    public static final int PROOF_RUNNING_TIME = 5*1000; // in ms

    public static void main(String[] args) {
        //solveInstance("data/TOP/gen3/m10/fnl4461_gen3_m10.txt", TopFormat.KUL,1000);
        //solveInstance("data/TOP/gen3/m3/d1291_gen3_m3.txt", TopFormat.KUL,1000);
        solveInstance("data/TOP/gen3/m3/u2319_gen3_m3.txt", TopFormat.KUL,1000);
        //solveInstance("data/TOP/Chaos/p4.2.b.txt", TopFormat.KUL,1000);
        //solveInstance("data/TOP/Chaos/p4.2.l.txt", TopFormat.KUL,1000);

        // solve every Chaos instance, then print a summary
        /*File[] files = new File("data/TOP/Chaos").listFiles();
        Arrays.sort(files);
        StringBuilder summary = new StringBuilder();
        for (File file : files) {
            try {
                summary.append(file.getName()).append(": ").append(solveInstance(file.getPath(), TopFormat.KUL, 1000)).append("\n");
            } catch (Exception e) {
                System.out.println("Exception for : "+ file.getPath() +  " : " + e.getMessage());
            }
        }

        System.out.println("===================== Summary =====================");
        System.out.print(summary);*/

    }

    public static String solveInstance(String instancePath, TopFormat format, int scaling) {
        TOPInstance instance = new TOPInstance(instancePath, format, scaling);
        System.out.println("running : " +  instancePath);

        // no vehicle can go from the start depot to the end depot: no route is possible, the best profit is 0
        if (instance.distanceMatrix[0][instance.nVert-1] > instance.maxTime) {
            return "0 (optimal)";
        }

        // ===================== Decision variables =====================
        CPSolver solver = makeSolver();
        CPSeqVar[] routes = new CPSeqVar[instance.nPath]; // All the routes (one for each vehicle)
        CPIntVar[] distance = makeIntVarArray(solver, instance.nPath, 0, instance.maxTime); // Store the total distance of a route.
        CPIntVar totalProfit = makeIntVar(solver, 0, instance.maxProfit);

        for (int i = 0; i < instance.nPath; i++) {
            routes[i] = makeSeqVar(solver, instance.nVert, 0, instance.nVert-1);
        }

        // ===================== Constrains =====================
        for(int i = 0; i < instance.nPath; i++) {
            solver.post(new Distance(routes[i], instance.distanceMatrix, distance[i], false)); // triangular inequality holds by construction
        }

        // A node can be visited once or none across all routes except for node 0 and node n-1
        for (int node = 1; node < instance.nVert-1; node++) {
            CPIntVar[] visits = new CPIntVar[instance.nPath];
            for (int route = 0; route < instance.nPath; route++) {
                visits[route] = routes[route].getNodeVar(node).isRequired();
            }
            solver.post(new Sum(visits, makeIntVar(solver, 0, 1)));
        }

        // Total profit: profit of each customer (depots excluded) visited by a route
        CPIntVar[] nodeProfit = new CPIntVar[instance.nPath * (instance.nVert - 2)];
        int k = 0;
        for (int i = 0; i < instance.nPath; i++) {
            for (int j = 1; j < instance.nVert - 1; j++) {
                nodeProfit[k++] = mul(routes[i].getNodeVar(j).isRequired(), instance.profit[j]);
            }
        }
        solver.post(new Sum(nodeProfit, totalProfit));

        CPIntVar totalDistance = sum(distance);
        //int m = Math.multiplyExact(instance.nPath, instance.maxTime) + 1;
        //Math.multiplyExact(instance.maxProfit, m); // throws if the objective overflows: lower the scaling
        //CPIntVar objective = sum(mul(totalProfit, m), minus(totalDistance));
        CPIntVar objective = totalProfit;

        // ===================== Search =====================
        DFSearch search = makeDfs(solver, firstFailBinary(routes));

        int[][] bestRoutes = new int[instance.nPath][instance.nVert];
        int[] bestLength = new int[instance.nPath];
        int[] bestProfit = new int[1];

        search.onSolution(() -> {
            for (int i = 0; i < instance.nPath; i++)
                bestLength[i] = routes[i].fillNode(bestRoutes[i], MEMBER_ORDERED);
            bestProfit[0] = totalProfit.min();
            System.out.println("Solution found: with total profit: " + totalProfit + " and total distance: " + totalDistance);
        });

        System.out.println("--- Running solver ---");
        Objective objectiveFunction = solver.maximize(objective);
        search.optimize(objectiveFunction, s -> s.numberOfSolutions() >= 1);

        // ===================== LNS =====================
        Random random = new Random(42);

        long start = System.currentTimeMillis();
        long end = start + LNS_RUNNING_TIME;

        int nCustomers = instance.nVert - 2;
        int nIter = 0;
        int nCompleted = 0;
        while (System.currentTimeMillis() < end) {
            int nCustomerDestroyed = random.nextInt(
                    Math.max(1, nCustomers * MIN_CUSTOMER_DESTROYED / 100),
                    Math.max(1, nCustomers * MAX_CUSTOMER_DESTROYED / 100) + 1);
            int nSequenceCustomerDestroyed = random.nextInt(
                    Math.max(1, nCustomers * MIN_SEQ_CUSTOMER_DESTROYED / 100),
                    Math.max(1, nCustomers * MAX_SEQ_CUSTOMER_DESTROYED / 100) + 1);

            int destroyId = random.nextInt(DESTROYS.length);
            Destroy destroy = DESTROYS[destroyId];

            boolean[] relaxed = switch (destroy) {
                case WORST_COST_SEQUENCE ->
                        worstCostSequenceDestroy(instance, bestRoutes, bestLength, nCustomerDestroyed, nSequenceCustomerDestroyed, random);
                case WORST_PROFIT_SEQUENCE ->
                        worstProfitSequenceDestroy(instance, bestRoutes, bestLength, nCustomerDestroyed, nSequenceCustomerDestroyed, random);
                case RANDOM_SEQUENCE ->
                        randomSequenceDestroy(instance, bestRoutes, bestLength, nCustomerDestroyed, nSequenceCustomerDestroyed, random);
                case SHAW ->
                        shawDestroy(instance, bestRoutes, bestLength, nCustomerDestroyed, nSequenceCustomerDestroyed, random);
                case RECTANGLE ->
                        rectangleDestroy(instance, bestRoutes, bestLength, nCustomerDestroyed, nSequenceCustomerDestroyed, random);
                default ->
                        worstCostEfficiencySequenceDestroy(instance, bestRoutes, bestLength, nCustomerDestroyed, nSequenceCustomerDestroyed, random);
            };

            SearchStatistics stats = search.optimizeSubjectTo(objectiveFunction, s -> s.numberOfFailures() >= 3_000, () -> {
                // keep the non-relaxed customers in their route and order
                for (int i = 0; i < instance.nPath; i++) {
                    int current = bestRoutes[i][0];
                    for (int j = 1; j < bestLength[i] - 1; j++) {
                        int node = bestRoutes[i][j];
                        if (!relaxed[node]) {
                            solver.post(insert(routes[i], current, node));
                            current = node;
                        }
                    }
                }
            });
            nIter++;
            if (stats.isCompleted())
                nCompleted++;
        }

        System.out.printf("LNS: %d neighborhoods, %d fully explored%n", nIter, nCompleted);

        // ===================== Optimal search =====================
        long proofEnd = System.currentTimeMillis() + PROOF_RUNNING_TIME;
        SearchStatistics proof = search.optimize(objectiveFunction, s -> System.currentTimeMillis() >= proofEnd);
        String result = bestProfit[0] + (proof.isCompleted() ? " (optimal)" : " (not proven)");
        System.out.println("best profit: " + result);
        return result;
    }

    // ===================== Destroy operators =====================
    public enum Destroy {
        WORST_COST_SEQUENCE,
        WORST_PROFIT_SEQUENCE,
        RANDOM_SEQUENCE,
        SHAW,
        RECTANGLE,
        WORST_COST_EFFICIENCY_SEQUENCE
    }

    private static final Destroy[] DESTROYS = Destroy.values();
    private static final double SELECTION_DETERMINISM = 3;

    public static boolean[] worstCostSequenceDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen, Random random) {
        return worstSequenceDestroy(instance, routes, length, q, maxSeqLen, random,
                (route, from, to) -> windowCost(instance, route, from, to));
    }

    public static boolean[] worstProfitSequenceDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen, Random random) {
        return worstSequenceDestroy(instance, routes, length, q, maxSeqLen, random,
                (route, from, to) -> -windowProfit(instance, route, from, to));
    }

    public static boolean[] randomSequenceDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen, Random random) {
        return worstSequenceDestroy(instance, routes, length, q, maxSeqLen, random,
                (route, from, to) -> random.nextDouble());
    }

    public static boolean[] worstCostEfficiencySequenceDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen, Random random) {
        return worstSequenceDestroy(instance, routes, length, q, maxSeqLen, random,
                (route, from, to) -> -costEfficiency(instance, route, from, to));
    }

    public static boolean[] shawDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen, Random random) {
        boolean[] relaxed = new boolean[instance.nVert];
        List<Integer> visited = visitedCustomers(routes, length);
        if (visited.isEmpty())
            return relaxed;
        List<Integer> removed = new ArrayList<>();
        removed.add(visited.remove(random.nextInt(visited.size())));
        while (removed.size() < q && !visited.isEmpty()) {
            int ref = removed.get(random.nextInt(removed.size()));
            visited.sort(Comparator.comparingInt(j -> instance.distanceMatrix[ref][j]));
            removed.add(visited.remove(select(random, visited.size())));
        }
        for (int node : removed)
            relaxed[node] = true;
        return relaxed;
    }

    public static boolean[] rectangleDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen, Random random) {
        boolean[] relaxed = new boolean[instance.nVert];
        List<Integer> visited = visitedCustomers(routes, length);
        int target = Math.min(q, visited.size());
        int seed = 1 + random.nextInt(instance.nVert - 2);
        int width = Arrays.stream(instance.xCoord).max().getAsInt() - Arrays.stream(instance.xCoord).min().getAsInt();
        int height = Arrays.stream(instance.yCoord).max().getAsInt() - Arrays.stream(instance.yCoord).min().getAsInt();
        long step = Math.max(1, Math.max(width, height) / 20);
        int nInside = 0;
        for (long halfSide = step; nInside < target; halfSide += step) {
            nInside = 0;
            for (int node : visited) {
                long dx = Math.abs(instance.xCoord[node] - instance.xCoord[seed]);
                long dy = Math.abs(instance.yCoord[node] - instance.yCoord[seed]);
                relaxed[node] = dx <= halfSide && dy <= halfSide;
                if (relaxed[node])
                    nInside++;
            }
        }
        return relaxed;
    }

    @FunctionalInterface
    private interface SequenceScore {
        double score(int[] route, int from, int to);
    }

    private static boolean[] worstSequenceDestroy(TOPInstance instance, int[][] routes, int[] length, int q, int maxSeqLen,
                                                  Random random, SequenceScore score) {
        record Sequence(int route, int from, int to, double score) {}
        boolean[] relaxed = new boolean[instance.nVert];
        int nRelaxed = 0;
        while (nRelaxed < q) {
            List<Sequence> candidates = new ArrayList<>();
            // fall back to shorter sequences if no sequence of the drawn length is left
            for (int seqLen = random.nextInt(1, Math.min(maxSeqLen, q - nRelaxed) + 1); seqLen >= 1 && candidates.isEmpty(); seqLen--) {
                for (int r = 0; r < routes.length; r++) {
                    for (int from = 1; from + seqLen - 1 <= length[r] - 2; from++) {
                        int to = from + seqLen - 1;
                        if (isIntact(relaxed, routes[r], from, to))
                            candidates.add(new Sequence(r, from, to, score.score(routes[r], from, to)));
                    }
                }
            }
            if (candidates.isEmpty())
                break; // every visited customer is already removed
            candidates.sort(Comparator.comparingDouble(Sequence::score).reversed());
            Sequence chosen = candidates.get(select(random, candidates.size()));
            for (int i = chosen.from(); i <= chosen.to(); i++)
                relaxed[routes[chosen.route()][i]] = true;
            nRelaxed += chosen.to() - chosen.from() + 1;
        }
        return relaxed;
    }

    private static boolean isIntact(boolean[] relaxed, int[] route, int from, int to) {
        for (int i = from; i <= to; i++)
            if (relaxed[route[i]])
                return false;
        return true;
    }

    private static int select(Random random, int size) {
        return (int) (Math.pow(random.nextDouble(), SELECTION_DETERMINISM) * size);
    }

    private static List<Integer> visitedCustomers(int[][] routes, int[] length) {
        List<Integer> visited = new ArrayList<>();
        for (int r = 0; r < routes.length; r++)
            for (int i = 1; i < length[r] - 1; i++)
                visited.add(routes[r][i]);
        return visited;
    }

    private static int windowCost(TOPInstance instance, int[] route, int from, int to) {
        int[][] t = instance.distanceMatrix;
        int cost = t[route[from - 1]][route[from]] + t[route[to]][route[to + 1]] - t[route[from - 1]][route[to + 1]];
        for (int i = from; i < to; i++)
            cost += t[route[i]][route[i + 1]];
        return cost;
    }

    private static int windowProfit(TOPInstance instance, int[] route, int from, int to) {
        int profit = 0;
        for (int i = from; i <= to; i++)
            profit += instance.profit[route[i]];
        return profit;
    }

    private static double costEfficiency(TOPInstance instance, int[] route, int from, int to) {
        int cost = windowCost(instance, route, from, to);
        return cost <= 0 ? Double.POSITIVE_INFINITY : (double) windowProfit(instance, route, from, to) / cost;
    }
}
