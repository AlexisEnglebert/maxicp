#!/bin/sh
# Experiments of the FDSSeqVar paper (Tables 3 and 4), about 15 minutes on a laptop.
# Usage (from the maxicp root): sh scripts/fdsseqvar_benchmarks.sh > results.txt
# Each line: instance search timeLimitSec nSol best tFirst tBest tEnd failures completed (times in ms).
set -e
mvn -q compile -Djacoco.skip=true
mvn -q dependency:build-classpath -Dmdep.outputFile=target/cp.txt
CP="target/classes:$(cat target/cp.txt)"
run() { java -cp "$CP" org.maxicp.cp.examples.modeling.FDSSeqVarBenchmark "$@" | tail -1; }

# Table 3: optimality proofs (safety limit 300 s); the second of two runs in the same JVM is kept (JIT warm-up)
for i in n40w60.001 n40w80.001 n40w100.001 pr01-7 pr01-8; do
  for s in nary primal uniform pool fds; do run $i $s 300 2; done
done

# Table 4: anytime runs on Cordeau's DARP instances
for i in pr01 pr02 pr03 pr04; do
  for s in nary uniform fds; do run $i $s 5; run $i $s 20; done
  run $i primal 5
  run $i bb 5
done
