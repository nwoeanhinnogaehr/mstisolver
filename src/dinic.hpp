#pragma once

#include <limits>
#include <vector>
using namespace std;

// This max flow code is modified from the kactl code archive
// at https://github.com/kth-competitive-programming/kactl/
// the original author information is below.

/**
 * Author: chilli
 * Date: 2019-04-26
 * License: CC0
 * Source: https://cp-algorithms.com/graph/dinic.html
 * Description: Flow algorithm with complexity $O(VE\log U)$ where $U = \max |\text{cap}|$.
 * $O(\min(E^{1/2}, V^{2/3})E)$ if $U = 1$; $O(\sqrt{V}E)$ for bipartite matching.
 * Status: Tested on SPOJ FASTFLOW and SPOJ MATCHING, stress-tested
 */
struct MaxFlow {
    struct DEdge {
        int to, rev;
        int64_t c, oc;
    };
    vector<int> lvl, ptr, q;
    vector<vector<DEdge>> adj;
	vector<int64_t> visit;
	int64_t iter = 1;
    MaxFlow(int n) : lvl(n), ptr(n), q(n), adj(n), visit(n) {}
    void addEdge(int a, int b, int64_t c, int64_t rcap = 0) {
        adj[a].push_back({ b, (int)adj[b].size(), c, c });
        adj[b].push_back({ a, (int)adj[a].size() - 1, rcap, rcap });
    }
	int64_t dfs_ff(int v, int t, int64_t f) {
		if (visit[v] == iter) return 0;
		if (v == t) return f;
		visit[v] = iter;
		for (DEdge &e : adj[v])
			if (e.c > 0)
				if (int p = dfs_ff(e.to, t, min(f, e.c))) {
					e.c -= p;
					adj[e.to][e.rev].c += p;
					return p;
				}
		return 0;
	}
	pair<int64_t, bool> calc_ff(int s, int t, int64_t limit=numeric_limits<int64_t>::max()) {
		int64_t flow = 0;
        bool cut = false;
		while (int p = dfs_ff(s, t, numeric_limits<int64_t>::max())) {
			iter++, flow += p;
            if (flow >= limit) {
                cut = true;
                break;
            }
        }
		iter++;
		return {flow, cut};
	}
    int64_t dfs(int v, int t, int64_t f) {
        if (v == t || !f)
            return f;
        for (int& i = ptr[v]; i < adj[v].size(); i++) {
            DEdge& e = adj[v][i];
            if (lvl[e.to] == lvl[v] + 1)
                if (int p = dfs(e.to, t, min(f, e.c))) {
                    e.c -= p, adj[e.to][e.rev].c += p;
                    return p;
                }
        }
        return 0;
    }
    int64_t calc(int s, int t, int64_t limit=numeric_limits<int64_t>::max()) {
        if (s == t) return numeric_limits<int64_t>::max();

        // do 1 iter of ford fulkerson; if that solves it then we're done early
        auto [ff, cut] = calc_ff(s, t, min(1l, limit));
        if (!cut || (cut && ff >= limit))
            return ff;

        int64_t flow = ff;
        q[0] = s;
        do {
            lvl.assign(q.size(), 0);
            ptr.assign(q.size(), 0);
            int qi = 0, qe = lvl[s] = 1;
            while (qi < qe && !lvl[t]) {
                int v = q[qi++];
                for (DEdge e : adj[v])
                    if (!lvl[e.to] && e.c)
                        q[qe++] = e.to, lvl[e.to] = lvl[v] + 1;
            }
            while (int p = dfs(s, t, numeric_limits<int64_t>::max())) {
                flow += p;
                if (flow >= limit)
                    return flow;
            }
        } while (lvl[t]);
        return flow;
    }
    void reset_flow() {
        for (int i = 0; i < adj.size(); i++)
            for (DEdge &e : adj[i])
                e.c = e.oc;
    }
};
