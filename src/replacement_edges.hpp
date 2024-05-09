#pragma once

#include "instance.hpp"
#include "union_find.hpp"
#include <vector>

struct ReplacementEdges {
    ReplacementEdges(Instance &inst, vector<bool> cur_x) : cur_x(cur_x), inst(inst) {}
    void compute() {
        height.resize(inst.n());
        parent_edge.assign(inst.n(), -1);
        parent.assign(inst.n(), -1);
        rep.assign(inst.m(), -1);
        compute_mst();
        dfs(0, 0, -1);
        UnionFind uf(inst.n());
        for (int i : non_mst_edges) {
            Edge &e = inst.edge(i);
            int x = uf.find(e.s);
            int y = uf.find(e.t);
            while (x != y) {
                if (height[x] < height[y]) swap(x, y);
                ASSERT(e.id > parent_edge[x]);
                rep[parent_edge[x]] = e.id;
                uf.join(parent[x], x);
                x = uf.find(parent[x]);
            }
        }
    }
    bool is_in_mst(Edge &e) {
        return rep[e.id] >= 0;
    }
    Edge &replacement(Edge &e) {
        ASSERT(rep[e.id] >= 0);
        return inst.edge(rep[e.id]);
    }

private:
    void compute_mst() {
        in_mst.assign(inst.m(), 0);
        UnionFind uf(inst.n());
        for (Edge &e : inst.edges()) {
            if (cur_x[e.id]) continue;
            if (uf.find(e.s) != uf.find(e.t)) {
                uf.join(e.s, e.t);
                in_mst[e.id] = true;
            } else
                non_mst_edges.push_back(e.id);
        }
    }
    void dfs(int i, int p, int pe) {
        height[i] = height[p]+1;
        parent_edge[i] = pe;
        parent[i] = p;
        for (Edge &e : inst.adj(i)) {
            if (!in_mst[e.id]) continue;
            if (e.s == i && e.t != p) dfs(e.t, i, e.id);
            else if (e.t == i && e.s != p) dfs(e.s, i, e.id);
        }
    }
    
    vector<bool> in_mst;
    vector<int> non_mst_edges;
    vector<int> height;
    vector<int> parent;
    vector<int> parent_edge;
    vector<int> rep;
    vector<bool> cur_x;
    Instance &inst;
};
