#pragma once

#include "replacement_edges.hpp"
#include "output.hpp"
#include <cmath>
#include <chrono>
using namespace std::chrono;

struct Greedy {
    Greedy(Instance &inst) : inst(inst) {}

    Soln compute() {
        auto start_time = high_resolution_clock::now();
        Soln s;
        s.x.resize(inst.m());
        s.y.assign(inst.m(), 0);
        s.weight = 0;
        UnionFind uf(inst.n());
        for (Edge &e : inst.edges()) {
            if (s.x[e.id]) continue;
            if (uf.find(e.s) != uf.find(e.t)) {
                s.y[e.id] = true;
                s.weight += e.lo_wt;
                uf.join(e.s, e.t);
            }
        }

        int64_t used = 0;
        while (true) {
            ReplacementEdges re(inst, s.x);
            re.compute();
            double best = -1;
            int best_edge = -1;
            for (int j = 0; j < inst.k(); j++) {
                Edge &e = inst.edge(j);
                if (used + e.up_wt > inst.cap())
                    continue;
                if (!s.y[e.id]) continue;

                int64_t gain = 0;
                int64_t cost = 0;
                Soln sol_chain = s;
                int64_t used_chain = used;
                Edge f = e;
                while (true) {
                    if (sol_chain.x[f.id])
                        break;
                    if (used_chain + f.up_wt > inst.cap())
                        break;
                    ReplacementEdges re_chain(inst, sol_chain.x);
                    re_chain.compute();
                    gain = re_chain.replacement(f).lo_wt - e.lo_wt;
                    cost += f.up_wt;
                    double val = gain / double(max(1l, cost));
                    if (val > best) {
                        best = val;
                        best_edge = e.id;
                    }
                    sol_chain.x[f.id] = true;
                    used_chain += f.up_wt;
                    f = re_chain.replacement(f);
                }
            }
            if (best_edge == -1)
                break;
            Edge &e = inst.edge(best_edge);
            used += e.up_wt;
            s.x[best_edge] = true;
            s.y[best_edge] = false;
            s.y[re.replacement(e).id] = true;
            s.weight += re.replacement(e).lo_wt - e.lo_wt;
        }
        auto end_time = high_resolution_clock::now();
        OUT("greedy_time " << duration<double, std::milli>(end_time - start_time).count() << std::endl);
        return s;
    }

private:
    Instance &inst;
};
