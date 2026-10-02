package org.maxicp.cp.examples.utils;

import org.maxicp.modeling.SeqVar;
import org.maxicp.modeling.algebra.sequence.SeqStatus;
import org.maxicp.util.algo.DistanceMatrix;
import org.maxicp.util.io.InputReader;
import java.util.Arrays;
import java.util.List;

import static org.maxicp.util.algo.DistanceMatrix.enforceTriangularInequality;

public class TOPInstance {

    public int nVert;
    public int nPath;
    public int maxTime;
    public int[][] distanceMatrix;
    public int[] profit;
    public int scaling;

    public int maxTotalDistances;
    public int maxProfit;

    public TOPInstance(String path, int scaling) {
        InputReader reader = new InputReader(path);
        this.scaling = scaling;

        reader.getString();
        nVert = reader.getInt();

        reader.getString();
        nPath = reader.getInt();

        reader.getString();
        maxTime = (int)(reader.getDouble() * this.scaling);

        maxTotalDistances = 0;
        maxProfit = 0;

        int[] xCoord = new int[nVert];
        int[] yCoord = new int[nVert];

        profit = new int[nVert];
        distanceMatrix = new int[nVert][nVert];

        // Assume that reader.getDouble() * this.scaling will not overflow.
        for(int i = 0; i < nVert; i++) {
            xCoord[i] = (int)(reader.getDouble() * this.scaling);
            yCoord[i] = (int)(reader.getDouble() * this.scaling);
            profit[i] = reader.getInt() * this.scaling;
        }

        for(int i = 0; i < nVert; i++) {
            for(int j = i+1; j < nVert; j++) {
                int dist = Math.round((float)Math.sqrt(Math.pow(xCoord[i]-xCoord[j], 2) + Math.pow(yCoord[i]-yCoord[j],2)));
                distanceMatrix[i][j] =  dist;
                distanceMatrix[j][i] =  dist;
            }
        }

        enforceTriangularInequality(distanceMatrix);

        for(int i = 0; i < nVert; i++) {
            maxTotalDistances += distanceMatrix[i][0];
            maxProfit += profit[i];
        }
    }

    public String toString() {
        return "Instance {" +
                "\nnVert = " + nVert +
                "\nnPath = " + nPath +
                "\nmaxTime = " + maxTime +
                "\ndistanceMatrix = " + Arrays.deepToString(distanceMatrix) +
                "\nprofit = " + Arrays.toString(profit) +
                "\n}";
    }

    /**
     * Checks that the given routes form a valid TOP solution:
     * - exactly nPath routes;
     * - each route goes from node 0 (start depot) to node nVert (end depot copy);
     * - depots never appear inside a route;
     * - each customer is visited at most once over all routes;
     * - the length of each route is at most maxTime.
     *
     * Must be called while the sequences are fixed (e.g. inside onSolution).
     *
     * @param result the routes, one per vehicle
     * @return true if the solution is valid
     */
    public boolean checkSolution(SeqVar[] result) {
        return checkSolution(result, -1);
    }

    /**
     * Same as {@link #checkSolution(SeqVar[])}, but also checks that the collected profit
     * equals {@code expectedProfit} (ignored if negative).
     */
    public boolean checkSolution(SeqVar[] result, int expectedProfit) {
        if (result == null || result.length != nPath) {
            System.err.printf("[CHECK] expected %d routes, got %d%n", nPath, result == null ? 0 : result.length);
            return false;
        }
        int startDepot = 0;
        int endDepot = nVert - 1;
        boolean[] visited = new boolean[nVert + 1];
        int[] nodes = new int[nVert + 1];
        int totalProfit = 0;

        for (int r = 0; r < nPath; r++) {
            SeqVar route = result[r];
            if (!route.isFixed()) {
                System.err.printf("[CHECK] route %d is not fixed%n", r);
                return false;
            }
            int n = route.fillNode(nodes, SeqStatus.MEMBER_ORDERED);
            if (n < 2 || nodes[0] != startDepot || nodes[n - 1] != endDepot) {
                System.err.printf("[CHECK] route %d must go from %d to %d: %s%n",
                        r, startDepot, endDepot, Arrays.toString(Arrays.copyOf(nodes, n)));
                return false;
            }
            int length = 0;
            for (int i = 0; i < n; i++) {
                int node = nodes[i];
                if (node < 0 || node > endDepot) {
                    System.err.printf("[CHECK] route %d contains invalid node %d%n", r, node);
                    return false;
                }
                if (i > 0) {
                    length += distanceMatrix[nodes[i - 1]][node];
                }
                if (i == 0 || i == n - 1)
                    continue; // depots, already checked
                if (node == startDepot || node == endDepot) {
                    System.err.printf("[CHECK] route %d visits a depot in the middle (position %d)%n", r, i);
                    return false;
                }
                if (visited[node]) {
                    System.err.printf("[CHECK] node %d is visited more than once (seen again in route %d)%n", node, r);
                    return false;
                }
                visited[node] = true;
                totalProfit += profit[node];
            }
            if (length > maxTime) {
                System.err.printf("[CHECK] route %d has length %d > maxTime %d%n", r, length, maxTime);
                return false;
            }
        }
        if (expectedProfit >= 0 && totalProfit != expectedProfit) {
            System.err.printf("[CHECK] collected profit %d differs from expected %d%n", totalProfit, expectedProfit);
            return false;
        }
        return true;
    }
}