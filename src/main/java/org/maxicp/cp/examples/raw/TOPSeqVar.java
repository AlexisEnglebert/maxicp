package org.maxicp.cp.examples.raw;

import org.maxicp.cp.engine.constraints.LessOrEqual;
import org.maxicp.cp.engine.constraints.MulCte;
import org.maxicp.cp.engine.constraints.Sum;
import org.maxicp.cp.engine.constraints.seqvar.Distance;
import org.maxicp.cp.engine.core.CPIntVar;
import org.maxicp.cp.engine.core.CPSeqVar;
import org.maxicp.cp.engine.core.CPSolver;
import org.maxicp.cp.examples.utils.TOPInstance;
import org.maxicp.search.DFSearch;
import org.maxicp.search.Objective;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import static org.maxicp.cp.CPFactory.*;
import static org.maxicp.modeling.algebra.sequence.SeqStatus.*;
import static org.maxicp.search.Searches.*;

public class TOPSeqVar {

    public static final int  MIN_CUSTOMER_DESTROYED = 1;
    public static final int MAX_CUSTOMER_DESTROYED = 20;

    public static final int MIN_SEQ_CUSTOMER_DESTROYED = 1;
    public static final int MAX_SEQ_CUSTOMER_DESTROYED = 20;

    public static final int LNS_RUNNING_TIME = 100*1000; // in ms

    public static void main(String[] args) {
        solveInstance("data/TOP/Chaos/p7.3.n.txt", 1000);
    }
    public static void solveInstance(String instancePath, int scaling) {
        TOPInstance instance = new TOPInstance(instancePath, scaling);

        // ===================== Decision variables =====================
        CPSolver solver = makeSolver();
        CPIntVar maxTime = makeIntVar(solver, instance.maxTime+1);
        CPSeqVar[] routes = new CPSeqVar[instance.nPath]; // All the routes (one for each vehicle)
        CPIntVar[] distance = new CPIntVar[instance.nPath]; // Store the total distance of a route.
        CPIntVar unit = makeIntVar(solver, 2);
        CPIntVar totalProfit = makeIntVar(solver, instance.maxProfit);
        CPIntVar[][] routeNodeProfit = new CPIntVar[instance.nPath][instance.nVert];
        CPIntVar[] perRouteProfit = new CPIntVar[instance.nPath];

        unit.fix(1);

        for(int i = 0; i < instance.nPath; i++) {
            routes[i] = makeSeqVar(solver, instance.nVert, 0, instance.nVert-1);
            distance[i] = makeIntVar(solver, instance.maxTotalDistances);
            System.out.println(routes[i]);
        }

        // ===================== Constrains =====================
        System.out.println(instance);
        maxTime.fix(instance.maxTime);

        // Each route must be <= than maxTime.
        for(int i = 0; i < instance.nPath; i++) {
            solver.post(new Distance(routes[i], instance.distanceMatrix, distance[i])); // Lui doit faire partis des troubleurs mais c'est le timide.
            solver.post(new LessOrEqual(distance[i], maxTime));
        }

        // A node can be visited once or none across all routes except for node 0 and node n-1
        for (int node = 1; node < instance.nVert-1; node++) {
            CPIntVar[] visits = new CPIntVar[instance.nPath];
            CPIntVar TotalVisit = makeIntVar(solver, instance.nVert);
            for (int route = 0; route < instance.nPath; route++) {
                visits[route] = routes[route].getNodeVar(node).isRequired();

            }
            solver.post(new Sum(visits, TotalVisit));
            solver.post(new LessOrEqual(TotalVisit, unit));
        }

        // Bind the profit to per route profit
        for (int i = 0; i < instance.nPath; i++) {
            for(int j = 0; j < instance.nVert; j++) {
                routeNodeProfit[i][j] = makeIntVar(solver, instance.maxProfit);
                solver.post(new MulCte(routes[i].getNodeVar(j).isRequired(), instance.profit[j], routeNodeProfit[i][j]));
            }
        }

        System.out.println(Arrays.deepToString(routeNodeProfit));
        for(int i = 0; i < instance.nPath; i++) {
            perRouteProfit[i] = makeIntVar(solver, instance.maxProfit);
            solver.post(new Sum(routeNodeProfit[i], perRouteProfit[i]));
        }
        solver.post(new Sum(perRouteProfit, totalProfit));

        // ===================== Search =====================

        DFSearch search = makeDfs(solver, firstFailBinary(routes));

        int[][] bestRoutes = new int[instance.nPath][instance.nVert];
        int[] bestLength = new int[instance.nPath];
        int[] bestProfit = new int[1];

        search.onSolution(() -> {
            for (int i = 0; i < instance.nPath; i++)
                bestLength[i] = routes[i].fillNode(bestRoutes[i], MEMBER_ORDERED);
            bestProfit[0] = totalProfit.min();
            System.out.println("yeaaah: " + Arrays.toString(routes) + " with total profit: " + totalProfit);
        });

        Objective objectiveFunction = solver.maximize(totalProfit);
        search.optimize(objectiveFunction, s -> s.numberOfSolutions() >= 1);

        // ===================== LNS =====================
        Random random = new Random(42);

        long start = System.currentTimeMillis();
        long end = start + LNS_RUNNING_TIME;

        while (System.currentTimeMillis() < end) {
            int nCustomerDestroyed = random.nextInt(MIN_CUSTOMER_DESTROYED, MAX_CUSTOMER_DESTROYED+1);
            int nSequenceCustomerDestroyed = random.nextInt(MIN_SEQ_CUSTOMER_DESTROYED, MAX_SEQ_CUSTOMER_DESTROYED+1);

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

            search.optimizeSubjectTo(objectiveFunction, s -> s.numberOfFailures() >= 4000, () -> {
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
        }

        System.out.println("best profit: " + bestProfit[0]);
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
