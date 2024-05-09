#pragma once

#include <set>

#include "instance.hpp"
#include "union_find.hpp"

struct InterdictionStructure {
    InterdictionStructure(Instance &inst)
        : uf(inst.n()), inst(inst) { }

    void init() {
        sol.x.resize(inst.m());
        sol.y.resize(inst.m());
        compute_mst();
    }

    void next() {
        ASSERT(mst_iter != mst.end());
        Edge e = inst.edge(*mst_iter);
        int a = uf.find(e.s);
        int b = uf.find(e.t);
        uf.joinBySize(a, b);
        mst_iter++;
    }
    void prev() {
        ASSERT(mst_iter != mst.begin());
        uf.undo();
        mst_iter--;
    }
    Edge &cur_edge() {
        return inst.edge(*mst_iter);
    }

    bool interdict() {
        int i = *mst_iter;
        ASSERT(sol.y[i]);
        ASSERT(!sol.x[i]);
        Edge e = inst.edge(i);

        int j = i+1;
        if (j == inst.m())
            return false;
        int undo = 0;
        for (; j < inst.m(); j++) {
            ASSERT(!sol.x[j]);
            Edge &f = inst.edge(j);
            int a = uf.find(f.s);
            int b = uf.find(f.t);
            if (sol.y[j]) {
                uf.joinBySize(a, b);
                undo++;
            } else if (a != b)
                break;
            if (j == inst.m()-1) {
                uf.undo(undo);
                return false;
            }
        }
        uf.undo(undo);

        sol.x[i] = true;
        sol.y[i] = false;
        sol.y[j] = true;
        sol.weight += inst.edge(j).lo_wt - e.lo_wt;
        mst.insert(j);
        mst_iter = mst.erase(mst_iter);
        stack.emplace_back(i, j);
        return true;
    }

    void un_interdict() {
        // undo the most recent interdiction
        auto [i, j] = stack.back();
        Edge &e = inst.edge(i);
        Edge &f = inst.edge(j);
        stack.pop_back();
        sol.x[i] = false;
        sol.y[i] = true;
        sol.y[j] = false;
        sol.weight += e.lo_wt - f.lo_wt;
        mst.erase(j);
        mst_iter = mst.insert(mst.end(), i);
    }

    Soln sol;

private:
    void compute_mst() {
        UnionFind uf(inst.n());
        sol.weight = 0;
        for (Edge &e : inst.edges()) {
            if (uf.find(e.s) != uf.find(e.t)) {
                uf.join(e.s, e.t);
                sol.weight += e.lo_wt;
                sol.y[e.id] = true;
                mst.insert(e.id);
            }
        }
        mst_iter = mst.begin();
    }
    set<int> mst;
    set<int>::iterator mst_iter;
    vector<pair<int,int>> stack;
    UnionFindWithUndo uf;
    
    // UnionFindWithUndo global_uf;
    Instance &inst;
};
