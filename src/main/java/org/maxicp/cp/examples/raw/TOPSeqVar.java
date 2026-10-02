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
import org.maxicp.search.SearchStatistics;

import java.util.Arrays;

import static org.maxicp.cp.CPFactory.*;
import static org.maxicp.search.Searches.*;

public class TOPSeqVar {

    public static void main(String[] args) {
        TOPInstance instance = new TOPInstance("data/TOP/Chaos/p4.2.b.txt", 1000);

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

        DFSearch search = makeDfs(solver, firstFailBinary(routes));

        search.onSolution(() -> {
            System.out.println("yeaaah: " + Arrays.toString(routes) + " with total profit: " + totalProfit);
            instance.checkSolution(routes);
        });

        Objective objectiveFunction = solver.maximize(totalProfit);
        SearchStatistics stats = search.optimize(objectiveFunction);
    }
}
