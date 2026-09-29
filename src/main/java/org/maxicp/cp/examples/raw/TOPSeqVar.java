package org.maxicp.cp.examples.raw;

import org.maxicp.cp.engine.constraints.Equal;
import org.maxicp.cp.engine.constraints.LessOrEqual;
import org.maxicp.cp.engine.constraints.Sum;
import org.maxicp.cp.engine.constraints.scheduling.CPSumCumulFunction;
import org.maxicp.cp.engine.constraints.seqvar.Distance;
import org.maxicp.cp.engine.core.CPIntVar;
import org.maxicp.cp.engine.core.CPSeqVar;
import org.maxicp.cp.engine.core.CPSolver;
import org.maxicp.cp.examples.utils.TOPInstance;
import org.maxicp.modeling.SeqVar;
import org.maxicp.search.DFSearch;
import org.maxicp.search.Searches;
import org.maxicp.util.algo.DistanceMatrix;

import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

import static org.maxicp.cp.CPFactory.*;
import static org.maxicp.modeling.algebra.sequence.SeqStatus.REQUIRED;
import static org.maxicp.search.Searches.*;

public class TOPSeqVar {

    public static void main(String[] args) {
        TOPInstance instance = new TOPInstance("data/TOP/Chaos/p4.2.a.txt", 1000);

        // ===================== Decision variables =====================
        CPSolver solver = makeSolver();
        CPIntVar maxTime = makeIntVar(solver, instance.maxTime+1);
        CPSeqVar[] routes = new CPSeqVar[instance.nPath]; // All the routes (one for each vehicle)
        CPIntVar[] distance = new CPIntVar[instance.nPath]; // Store the total distance of a route.

        for(int i = 0; i < instance.nPath; i++) {
            routes[i] = makeSeqVar(solver, instance.nVert+1, 0, instance.nVert);
            distance[i] = makeIntVar(solver, instance.maxTotalDistances);
            System.out.println(routes[i]);
        }

        // ===================== Constrains =====================
        maxTime.fix(instance.maxTime);

        // Each route must be <= than maxTime.
        for(int i = 0; i < instance.nPath; i++) {
            solver.post(new Distance(routes[i], instance.distanceMatrix, distance[i]));
            solver.post(new LessOrEqual(distance[i], maxTime));
        }

        // A node can be visited once across all routes except for node 0 and node n
        for (int node = 1; node < instance.nVert; node++) {
            CPIntVar[] visits = new CPIntVar[instance.nPath];
            for (int route = 0; route < instance.nPath; route++) {
                visits[route] = routes[route].getNodeVar(node).isRequired();
            }
            solver.post(new Sum(visits, 1));
        }

        DFSearch search = makeDfs(solver, firstFailBinary(routes));

        search.onSolution(() -> {
            System.out.println("yeaaah: " + Arrays.toString(routes));
        });
    }
}
