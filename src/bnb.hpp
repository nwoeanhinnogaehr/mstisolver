#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "instance.hpp"
#include "knapsack_ub.hpp"
#include "interdiction_structure.hpp"

struct BranchAndBound {
public:
    BranchAndBound(Instance &inst, Soln lb, KnapsackUpperBound *kp_ub, int64_t ub_ub = numeric_limits<int64_t>::min())
        : inst(inst), knapsack_ub(kp_ub), is(inst), lb(lb), ub_ub(ub_ub) {
    }
    Soln solve() {
        min_postfix.resize(inst.m());
        min_postfix[inst.m()-1] = inst.edge(inst.m()-1).up_wt;
        for (int i = inst.m()-2; i >= 0; i--)
            min_postfix[i] = min(min_postfix[i+1], inst.edge(i).up_wt);
        is.init();
        branch(0, inst.cap(), 0);
        return lb;
    }

private:
    void branch(int i, int64_t rem_cap, int prefix);
    InterdictionStructure is;
    Soln lb;
    vector<Soln> stack;
    KnapsackUpperBound *knapsack_ub;
    Instance &inst;
    int64_t ub_ub;
    uint64_t n_new_nodes = 0, n_new_pruned = 0;
    vector<int64_t> min_postfix;
};

void BranchAndBound::branch(int i, int64_t rem_cap, int prefix) {
    n_nodes++;
    n_new_nodes++;

    if (n_new_nodes == 1000000) {
        LOG("+ searched " << n_nodes/1000000 << "M nodes so far, pruned "
            << n_new_pruned << "/" << n_new_nodes << endl);
        LOG("  current node: ";
        for (int j = 0; j < min(is.cur_edge().id, 60); j++)
            cerr << is.sol.x[j];
        cerr << (is.cur_edge().id > 60 ? "..." : "") << endl);
        n_new_nodes = n_new_pruned = 0;
    }

    if (is.sol.weight > lb.weight) {
        LOG("+ improved lower bound to " << is.sol.weight << endl);
        lb = is.sol;
    }
    if (i == inst.n()-1 || is.cur_edge().id >= inst.k()) {
        return;
    }
    if (min_postfix[is.cur_edge().id] > rem_cap) {
        return;
    }

    if (knapsack_ub) {
        int64_t ub = knapsack_ub->get_bound(is.cur_edge().id, rem_cap, prefix);
        if (ub != numeric_limits<int64_t>::max() && ub + is.sol.weight <= max(lb.weight, ub_ub)) {
            n_pruned++;
            n_new_pruned++;
            return;
        }
    }

    Edge &e = is.cur_edge();
    ASSERT(is.sol.y[e.id]);
    if (e.up_wt <= rem_cap) {
        is.interdict();
        branch(i, rem_cap - e.up_wt, knapsack_ub && e.id < knapsack_ub->max_prefix_bits() ? prefix | (1<<e.id) : prefix);
        is.un_interdict();
    }
    is.next();
    branch(i+1, rem_cap, prefix);
    is.prev();
}
