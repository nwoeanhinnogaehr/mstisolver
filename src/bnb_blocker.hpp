#pragma once

#include <cstdint>
#include <limits>
#include <vector>
#include <unordered_map>

#include "instance.hpp"
#include "knapsack_ub.hpp"
#include "interdiction_structure.hpp"

struct BlockerBranchAndBound {
public:
    BlockerBranchAndBound(Instance &inst, KnapsackUpperBound *kp_ub)
        : inst(inst), knapsack_ub(kp_ub), is(inst) {
    }
    Soln solve() {
        incumbent = Soln(inst);
        best_cost = inst.min_cut() - 1;
        is.init();
        branch(0, 0, 0);
        return incumbent;
    }

private:
    void branch(int i, int64_t cur_cost, int prefix);
    InterdictionStructure is;
    Soln incumbent;
    KnapsackUpperBound *knapsack_ub;
    Instance &inst;
    int64_t best_cost;
    uint64_t n_new_nodes = 0, n_new_pruned = 0;
};

void BlockerBranchAndBound::branch(int i, int64_t cur_cost, int prefix) {
    n_nodes++; n_new_nodes++;

    if (n_new_nodes == 1000000) {
        LOG("+ searched " << n_nodes/1000000 << "M nodes so far, pruned "
            << n_new_pruned << "/" << n_new_nodes << endl);
        LOG("  current node: ";
        for (int j = 0; j < min(is.cur_edge().id, 60); j++)
            cerr << is.sol.x[j];
        cerr << (is.cur_edge().id > 60 ? "..." : "") << endl);
        n_new_nodes = n_new_pruned = 0;
    }

    if (is.sol.weight >= inst.r() && cur_cost <= best_cost) {
        LOG("+ improved incumbent to " << cur_cost << endl);
        best_cost = cur_cost;
        incumbent = is.sol;
    }
    if (i == inst.n()-1 || is.cur_edge().id >= inst.k() || is.sol.weight >= inst.r() || cur_cost > best_cost) {
        return;
    }

    Edge &e = is.cur_edge();
    ASSERT(is.sol.y[e.id]);
    if (knapsack_ub) {
        int64_t diff = best_cost - cur_cost;
        ASSERT(diff >= 0);
        ASSERT(diff <= inst.min_cut()-1);
        int64_t ub = knapsack_ub->get_bound(e.id, diff, prefix);
        if (ub != numeric_limits<int64_t>::max() && ub < inst.r() - is.sol.weight) {
            n_pruned++; n_new_pruned++;
            return;
        }
    }
    if (is.interdict()) {
        branch(i, cur_cost + e.up_wt, knapsack_ub && e.id < knapsack_ub->max_prefix_bits() ? prefix | (1<<e.id) : prefix);
        is.un_interdict();
    }
    is.next();
    branch(i+1, cur_cost, prefix);
    is.prev();
}
