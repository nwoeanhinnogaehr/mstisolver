#pragma once

#include <limits>
#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include <optional>
using namespace std;
using namespace std::chrono;

#include "assert.hpp"
#include "dinic.hpp"
#include "union_find.hpp"
#include "output.hpp"

struct Edge {
    Edge(int s, int t, int64_t lo_wt, int64_t up_wt)
        : s(s)
        , t(t)
        , lo_wt(lo_wt)
        , up_wt(up_wt)
        , id(-1) {}
    int s, t;
    int64_t lo_wt, up_wt;
    int id;
};
enum class ProblemType {
    MSTInterdiction = 0,
    EdgeBlocker
};
struct Instance {
public:
    Instance() {}
    // for mst interdiction problem
    Instance(int n, int m, int64_t cap, vector<Edge> edges, int k = -1, optional<int64_t> mincut = {}) {
        _n = n;
        _m = m;
        _cap = cap;
        _edges = edges;
        if (k == -1)
            _k = m;
        else
            _k = k;
        _min_cut = mincut;
        ASSERT(is_sorted(_edges.begin(), _edges.end(), [](const Edge& a, const Edge& b){
            return a.lo_wt < b.lo_wt;
        }));
        for (int i = 0; i < m; i++) {
            _edges[i].id = i;
        }
        _adj.assign(n, vector<Edge>());
        for (Edge& e : _edges) {
            _adj[e.s].push_back(e);
            _adj[e.t].push_back(e);
        }
    }
    // for edge blocker problem
    Instance(int n, int m, vector<Edge> edges) : Instance(n, m, -1, edges) { 
        _problem_type = ProblemType::EdgeBlocker;
        _cap = min_cut() - 1;
    }

    int n() { return _n; }
    int m() { return _m; }
    int64_t cap() {
        return _cap;
    }
    int k() {
        return _k;
    }
    int64_t r() {
        // ASSERT(_problem_type == ProblemType::EdgeBlocker);
        return _r;
    }
    void set_cap(int64_t cap) {
        _cap = cap;
    }
    void set_r(int64_t r) {
        _r = r;
    }
    ProblemType problem_type() { return _problem_type; }
    Edge &edge(int id) {
        return _edges[id];
    }
    vector<Edge> &edges() {
        return _edges;
    }
    vector<Edge> &adj(int v) {
        return _adj[v];
    }

    bool is_valid() {
        if (_problem_type == ProblemType::EdgeBlocker) {
            ASSERT(false);
            return false;
        } else {
            return min_cut() > cap();
        }
    }

    int64_t min_cut() {
        if (_min_cut) {
            return *_min_cut;
        } else {
            MaxFlow d(n());
            UnionFind uf(n());
            for (int i = 0; i < _m; i++) {
                Edge &e = edge(i);
                if (e.up_wt > 0) {
                    d.addEdge(e.s, e.t, e.up_wt, e.up_wt);
                    uf.join(e.s, e.t);
                }
            }
            for (int j = 1; j < n(); j++) {
                if (uf.find(0) != uf.find(j))
                    return _min_cut.emplace(0); // graph is disconnected
            }
            int64_t mincut = numeric_limits<int64_t>::max();
            for (int j = 1; j < n(); j++) {
                mincut = min(mincut, d.calc(0, j, mincut));
                d.reset_flow();
                if (mincut <= 1) break;
            }
            return _min_cut.emplace(mincut);
        }
    }

    Instance canonicalize() {
        // as described by Zenklusen, and later Linhares and Swamy
        auto start_time = high_resolution_clock::now();

        vector<int> breakpoints;
        breakpoints.push_back(0);
        for (int k = 1; k < m(); k++) {
            Edge &e = edge(k);
            if (e.lo_wt == edge(breakpoints.back()).lo_wt)
                continue;
            breakpoints.push_back(k);
        }
        breakpoints.push_back(m());

        int lo = 0;
        int hi = breakpoints.size()-1;
        while (lo < hi) {
            int mid = (lo + hi)/2;
            int breakpoint = breakpoints[mid];

            MaxFlow d(n());
            for (int j = 0; j < breakpoint; j++) {
                Edge &e = edge(j);
                d.addEdge(e.s, e.t, e.up_wt, e.up_wt);
            }
            bool cut = false;
            for (int j = 1; j < n(); j++) {
                if (d.calc(0, j, cap()+1) <= cap()) {
                    cut = true;
                    break;
                }
                d.reset_flow();
            }
            if (cut)
                lo = mid + 1;
            else
                hi = mid;
        }
        int K = breakpoints[lo]-1;
        vector<Edge> edges;
        for (int i = 0; i <= K; i++)
            edges.push_back(edge(i));
        ASSERT(lo >= 1);
        Instance c(n(), edges.size(), cap(), edges, breakpoints[lo-1], min_cut());

        auto end_time = high_resolution_clock::now();
        OUT("canonicalize_time " << duration<double, std::milli>(end_time - start_time).count() << std::endl);
        return c;
    }

    void write_to(ostream &os) {
        if (_problem_type == ProblemType::MSTInterdiction) {
            os << "problem_type msti" << endl;
            if (_cap >= 0)
                os << "cap " << _cap << endl;
        } else {
            os << "problem_type mebsp" << endl;
            if (_r >= 0)
                os << "target_weight" << _r << endl;
        }
        os << "n_verts " << _n << endl;
        os << "n_edges " << _m << endl;
        for (Edge &e : _edges)
            os << "edge " << e.s << " " << e.t << " " << e.lo_wt << " " << e.up_wt << endl;
    }
    static Instance read_from(istream &is) {
        int n = -1, m = -1;
        int64_t cap = -1, r = -1;
        ProblemType pt = ProblemType::MSTInterdiction;
        vector<Edge> edges;
        while (!is.eof()) {
            string s;
            is >> s;
            if (s == "n_verts") {
                is >> n;
            } else if (s == "n_edges") {
                is >> m;
            } else if (s == "cap") {
                is >> cap;
            } else if (s == "target_weight") {
                is >> r;
            } else if (s == "problem_type") {
                string ty;
                is >> ty;
                if (ty == "msti")
                    pt = ProblemType::MSTInterdiction;
                else if (ty == "mebsp")
                    pt = ProblemType::EdgeBlocker;
                else ASSERT(false);
            } else if (s == "edge") {
                int s, t;
                int64_t lo_wt, up_wt;
                is >> s >> t >> lo_wt >> up_wt;
                ASSERT(s >= 0 && t >= 0 && s < n && t < n);
                edges.push_back(Edge(s, t, lo_wt, up_wt));
            }
        }
        ASSERT(edges.size() == m);
        if (pt == ProblemType::MSTInterdiction) {
            ASSERT(cap >= 0);
            ASSERT(r == -1);
            return Instance(n, m, cap, edges);
        } else {
            ASSERT(cap == -1);
            Instance inst(n, m, edges);
            inst._r = r;
            return inst;
        }
    }

    // Generate a random instance with n vertices, m edges, lo_wt (up_wt)
    // selected uniformly at random in [1, max_lo_wt (max_up_wt)], and capacity
    // set to cap_scale times the max capacity that yields a valid instance
    // (with no leader cut).
    static Instance generate_random(int n, int m, int max_lo_wt, int max_up_wt,
                                    double cap_scale, unsigned int seed=-1) {
        std::random_device rd;
        std::mt19937 rng(seed == -1 ? rd() : seed);
        std::uniform_int_distribution<int> lo_wt_gen(1, max_lo_wt);
        std::uniform_int_distribution<int> up_wt_gen(1, max_up_wt);
        vector<Edge> edges;
        ASSERT(n >= 0 && m >= 0);
        ASSERT(m <= n * (n - 1) / 2);
        ASSERT(0.0 <= cap_scale && cap_scale <= 1.0);
        vector<Edge> possible_edges;
        // TODO: there are definitely faster ways to generate sparse graphs!
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                possible_edges.push_back(
                    Edge(i, j, lo_wt_gen(rng), up_wt_gen(rng)));
        shuffle(possible_edges.begin(), possible_edges.end(), rng);
        copy(possible_edges.begin(), possible_edges.begin()+m, back_inserter(edges));
        MaxFlow d(n);
        for (Edge &e : edges)
            d.addEdge(e.s, e.t, e.up_wt, e.up_wt);
        int64_t mincut = numeric_limits<int64_t>::max();
        for (int i = 1; i < n; i++) {
            mincut = min(mincut, d.calc(0, i, mincut));
            d.reset_flow();
        }
        int64_t cap = cap_scale * (mincut - 1);
        ASSERT(cap >= 0);
        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b){
            return a.lo_wt < b.lo_wt;
        });
        return Instance(n, m, cap, edges);
    }

    static Instance swamy_lowerbound(int n) {
        vector<Edge> edges;
        for (int i = 0; i < n-1; i++)
            edges.push_back(Edge(i, (i+1)%(n-1), 0, n));
        for (int i = 0; i < n-1; i++)
            edges.push_back(Edge(i, n-1, 1, 2*n));
        int cap = 2*n-2;
        return Instance(n, edges.size(), cap, edges);
    }

    static Instance long_comb_instance(int k, int alpha) {
        vector<Edge> edges;
        int inf = 2*k;
        for (int i = 0; i < k; i++) {
            int offset = 2*(k+2)*i;
            for (int j = 0; j <= k; j++) {
                edges.push_back(Edge(offset+j, offset+j+1, 0, inf));
                edges.push_back(Edge(offset+j+k+2, offset+j+k+3, 0, inf));
                if (j == k)
                    edges.push_back(Edge(offset+j+1, offset+j+k+3, alpha, inf));
                else
                    edges.push_back(Edge(offset+j+1, offset+j+k+3, j+1, 1));
            }
            if (i != k-1)
                edges.push_back(Edge(offset+k+2, offset+2*k+4, 0, inf));
        }
        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b){
            return a.lo_wt < b.lo_wt;
        });
        return Instance(2*(k+2)*k, edges.size(), 2*k-1, edges);
    }
    static Instance comb_instance(int k, int alpha) {
        vector<Edge> edges;
        int inf = 2*k;
        for (int i = 0; i < k; i++) {
            int offset = 6*i;
            edges.push_back(Edge(offset, offset+1, 0, inf));
            edges.push_back(Edge(offset+1, offset+2, 0, inf));
            edges.push_back(Edge(offset+3, offset+4, 0, inf));
            edges.push_back(Edge(offset+4, offset+5, 0, inf));
            edges.push_back(Edge(offset, offset+3, 0, k-1));
            edges.push_back(Edge(offset+1, offset+4, 1, 1));
            edges.push_back(Edge(offset+2, offset+5, alpha, inf));
            if (i != k-1)
                edges.push_back(Edge(offset+3, offset+6, 0, inf));
        }
        // edges.push_back(Edge(0, 6*(k-1)+3, 0, 1));
        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b){
            return a.lo_wt < b.lo_wt;
        });
        return Instance(6*k, edges.size(), 2*k-1, edges);
    }

protected:
    int _n, _m, _k;
    int64_t _cap, _r = -1;
    vector<Edge> _edges;
    vector<vector<Edge>> _adj;
    ProblemType _problem_type = ProblemType::MSTInterdiction;
    optional<int64_t> _min_cut;
};

struct Soln {
    vector<bool> x, y;
    int64_t weight;

    Soln() {}
    Soln(Instance &inst) : x(inst.m()), y(inst.m()), weight(0) {
        UnionFind uf(inst.n());
        for (Edge &e : inst.edges()) {
            if (uf.find(e.s) != uf.find(e.t)) {
                y[e.id] = 1;
                weight += e.lo_wt;
                uf.join(e.s, e.t);
            }
        }
    }
};
