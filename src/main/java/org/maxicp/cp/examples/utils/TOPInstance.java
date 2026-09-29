package org.maxicp.cp.examples.utils;

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

    public TOPInstance(String path, int scaling) {
        InputReader reader = new InputReader(path);
        this.scaling = scaling;

        reader.getString();
        nVert = reader.getInt();

        reader.getString();
        nPath = reader.getInt();

        reader.getString();
        maxTime = (int)(reader.getDouble() * this.scaling);


        int[] xCoord = new int[nVert];
        int[] yCoord = new int[nVert];

        profit = new int[nVert];
        int[][] distance = new int[nVert][nVert];

        // Assume that reader.getDouble() * this.scaling will not overflow.
        for(int i = 0; i < nVert; i++) {
            xCoord[i] = (int)(reader.getDouble() * this.scaling);
            yCoord[i] = (int)(reader.getDouble() * this.scaling);
            profit[i] = reader.getInt() * this.scaling;
        }

        for(int i = 0; i < nVert; i++) {
            for(int j = i+1; j < nVert; j++) {
                int dist = Math.round((float)Math.sqrt(Math.pow(xCoord[i]-xCoord[j], 2) + Math.pow(yCoord[i]-yCoord[j],2)) * this.scaling);
                distance[i][j] =  dist;
                distance[j][i] =  dist;
            }
        }

        // Duplicate the row 0 because seq var needs a distinct node for begin and end. So we add one row so that node 0 and n is the stard and end.
        distanceMatrix = DistanceMatrix.extendMatrixAtEnd(distance, List.of(0));
        enforceTriangularInequality(distanceMatrix);

        for(int i = 0; i < nVert; i++) {
            maxTotalDistances += distanceMatrix[i][0];
        }
    }

    public String toString() {
        return "Instance {" +
                "\nnVert = " + nVert +
                "\nnPath = " + nPath +
                "\nmaxTime = " + maxTime +
                "\ndistanceMatrix = " + Arrays.deepToString(distanceMatrix) +
                "\nproffit = " + Arrays.toString(profit) +
                "\n}";
    }
}
