#pragma once

#include "instance.hpp"
#include "dinic.hpp"
#include "assert.hpp"
#include "output.hpp"

#include <iostream>
#include <limits>
#include <vector>
#include <thread>
#include <atomic>
using namespace std;
using namespace std::chrono;

struct KnapsackUpperBound {
    virtual ~KnapsackUpperBound() {};
    virtual void compute_in_background() = 0;
    virtual int64_t get_bound(int i, int r, int prefix) = 0;
    virtual int baseline_weight() = 0;
    virtual int prefix_bits_computed() = 0;
    virtual int max_prefix_bits() = 0;
};

template <typename T>
struct KnapsackUpperBoundT : KnapsackUpperBound {
public:
    KnapsackUpperBoundT(Instance linst, int64_t scale, int max_bits = 30)
      : max_bits(max_bits)
      , scale(scale)
    {
        _dp_table.resize(max_bits);
        _edge_deltas.resize(max_bits);
        inst = linst;
        if (inst.problem_type() == ProblemType::MSTInterdiction) {
            cap = (inst.cap()+scale-1)/scale;
        } else {
            cap = (inst.min_cut()-1+scale-1)/scale;
        }
        _k = inst.k();
        lo_wt.resize(inst.m());
        up_wt.resize(inst.m());
        for (int i = 0; i < inst.m(); i++) {
            Edge &e = inst.edge(i);
            lo_wt[i] = e.lo_wt;
            up_wt[i] = e.up_wt/scale;
        }
        baseline = Soln(inst);
    }
    virtual ~KnapsackUpperBoundT() {
        _should_quit.store(true);
        _thread.join();
    }
    void compute_in_background() override {
        _thread = thread(&KnapsackUpperBoundT::compute_thread, this);
    }
    int64_t compute(int bits = 0) {
        ASSERT(_bits_computed.load() == bits-1);
        compute_edge_deltas(bits);
        _dp_table[bits].assign(uint64_t(cap+1)*_k*(1ull<<bits), numeric_limits<T>::max());

        for (int i = 0; i <= cap; i++)
            dp_rec(0, i, 0, bits);
        int64_t wt = dp_rec(0, cap, 0, bits) + baseline.weight;

        // The tables are not atomic to avoid unnecessary overhead,
        // so we need a fence here to ensure that other threads will
        // definitely see the updated tables when they read _bits_computed.
        atomic_thread_fence(memory_order_seq_cst);
        _bits_computed.store(bits);

        return wt;
    }
    bool is_computed() {
        return _bits_computed.load() >= 0;
    }
    int64_t get_bound(int i, int r, int prefix) override {
        int bits = prefix_bits_computed();
        if (bits == -1)
            return numeric_limits<int64_t>::max();
        if (bits > 0) {
            _dp_table[bits - 1].clear();
            _dp_table[bits - 1].shrink_to_fit();
            _edge_deltas[bits - 1].clear();
            _edge_deltas[bits - 1].shrink_to_fit();
        }
        int pre = prefix & ((1<<min(i, bits))-1);
        r = (r+scale-1)/scale;
        ASSERT(0 <= r && r <= cap);
        T entry = table(i, r, pre, bits);
        ASSERT(entry != numeric_limits<T>::max());
        return entry;
    }
    int baseline_weight() override {
        return baseline.weight;
    }
    int prefix_bits_computed() override {
        return _bits_computed.load(memory_order_seq_cst);
    }
    int max_prefix_bits() override {
        return max_bits;
    }

private:
    void compute_thread() {
        int64_t root_bound = numeric_limits<int64_t>::max();
        for (int bits = 0; bits < max_bits; bits++) {
            uint64_t size = uint64_t(cap+1)*_k*(1ull<<bits)*sizeof(T);
            LOG("+ computing dp upper bound with " << bits << " prefix bits..." << endl);
            if (3*size > MEM_LIMIT) { // 3 because it's the previous + next size of dp+delta tables, and previous is half the size
                LOG("+ hit memory limit for DP tables" << endl);
                return;
            }
            auto start_time = high_resolution_clock::now();
            root_bound = min(root_bound, compute(bits));
            if (_should_quit.load(memory_order_relaxed))
                return;
            OUT("dp_upper_bound_" << bits << " " << root_bound << endl);
            auto end_time = high_resolution_clock::now();
            OUT("dp_upper_bound_time_" << bits << " " << duration<double, std::milli>(end_time - start_time).count() << std::endl);
        }
        LOG("+ upper bound thread finished" << endl);
    }
    
    T dp_rec(int i, int r, int pre, int bits) {
        if (i >= _k)
            return 0;

        ASSERT(i >= 0);
        ASSERT(i < inst.m());
        ASSERT(r >= 0);
        ASSERT(r <= cap);

        Edge &e = inst.edge(i);
        T &entry = table(i, r, pre, bits);
        if (entry != numeric_limits<T>::max())
            return entry;

        if (_should_quit.load(memory_order_relaxed))
            return 0;

        if (up_wt[i] <= r) {
            if (i <= bits-1) { // in prefix
                    int upd_pre = pre | (1<<i);
                    entry = max(T(dp_rec(i+1, r-up_wt[i], upd_pre, bits) + edge_delta(i, r, pre, bits)),
                                dp_rec(i+1, r, pre, bits));
            } else // out of prefix
                    entry = max(T(dp_rec(i+1, r-up_wt[i], pre, bits) + edge_delta(i, r, pre, bits)),
                                dp_rec(i+1, r, pre, bits));
        } else
            entry = dp_rec(i+1, r, pre, bits);

        return entry;
    }

    void compute_edge_deltas(int bits) {
        _edge_deltas[bits].resize(uint64_t(1ull<<bits)*_k*(cap+1));

        // omp_set_num_threads(4);
        // #pragma omp parallel for
        for (int i = 0; i < _k; i++) {
            for (int pre = 0; pre < (1<<min(i, bits)); pre++) {
                Edge &e = inst.edge(i);
                int k = inst.m()-1;
                UnionFind uf(inst.n()), uf2(inst.n());
                for (int j = 0; j < bits; j++) {
                    Edge &f = inst.edge(j);
                    if (!(pre & (1<<j)) && uf.find(e.s) != uf.find(e.t)) {
                        uf.join(e.s, e.t);
                        uf2.join(e.s, e.t);
                    }
                }
                uf.join(e.s, e.t);
                for (int j = i+1; j < inst.m(); j++) {
                    Edge &f = inst.edge(j);
                    if (uf2.find(f.s) != uf2.find(f.t)) {
                        uf2.join(f.s, f.t);
                        if (uf.find(f.s) == uf.find(f.t)) {
                            k = j;
                            break;
                        }
                        uf.join(f.s, f.t);
                    }
                }
                for (int r = 0; r <= cap; r++)
                    edge_delta(e.id, r, pre, bits) = lo_wt[k]-lo_wt[i];
                int subtract_cost = 0;
                MaxFlow d(inst.n());
                for (int k = 0; k < e.id; k++) {
                    Edge &f = inst.edge(k);
                    if (k < bits) {
                        if (!(pre & (1<<k)))
                            d.addEdge(f.s, f.t, cap+1, cap+1);
                        else 
                            subtract_cost += up_wt[k];
                    } else
                        d.addEdge(f.s, f.t, up_wt[k], up_wt[k]);
                }
                ASSERT(e.s != e.t);
                int mincut = d.calc(e.s, e.t, cap+1);
                for (int r = max(0, cap - subtract_cost - mincut + 1); r <= cap; r++)
                    edge_delta(e.id, r, pre, bits) = 0;
                for (int j = e.id+1; j < inst.m(); j++) {
                    if (_should_quit.load(memory_order_relaxed))
                        return;
                    Edge e2 = inst.edge(j);
                    d.addEdge(e2.s, e2.t, cap+1, cap+1);
                    mincut += d.calc(e.s, e.t, cap+1-mincut);
                    for (int r = max(0, cap - subtract_cost - mincut + 1); r <= cap; r++)
                        edge_delta(e.id, r, pre, bits) = min(
                            edge_delta(e.id, r, pre, bits), T(lo_wt[j] - lo_wt[i]));
                    if (mincut > cap)
                        break;
                }
            }
        }
    }

    T &edge_delta(int i, int r, int pre, int bits) {
        ASSERT(pre <= ((1<<min(i,bits))-1));
        ASSERT(0 <= bits && bits < max_bits);
        return _edge_deltas[bits][uint64_t(pre) * (cap + 1) +
                                 uint64_t(i) * (cap + 1) * (1ull<<bits) + uint64_t(r)];
    }
    T &table(int i, int r, int pre, int bits) {
        ASSERT(pre <= ((1<<min(i,bits))-1));
        ASSERT(0 <= bits && bits < max_bits);
        return _dp_table[bits][uint64_t(pre) * (cap + 1) +
                              uint64_t(i * (cap + 1) * (1ull<<bits)) + uint64_t(r)];
    }
    Instance inst;
    Soln baseline;
    int cap;
    int _k;
    vector<T> lo_wt;
    vector<int64_t> up_wt;

    vector<vector<T>> _edge_deltas;
    vector<vector<T>> _dp_table;

    thread _thread;
    atomic<int> _bits_computed = -1;
    atomic<bool> _should_quit = false;

    int64_t scale = 1;
    int max_bits;
    uint64_t MEM_LIMIT = (1ull << 33) + (1ull << 32) + (1ull << 31); // 14GB limit
};

KnapsackUpperBound *make_knapsack_ub(Instance inst, int64_t scale, int max_bits) {
    uint64_t ub = inst.edge(inst.m()-1).lo_wt * (inst.n()-1); // weak
    if (ub < numeric_limits<uint8_t>::max()) {
        LOG("+ using 8bit DP table" << endl);
        return new KnapsackUpperBoundT<uint8_t>(inst, scale, max_bits);
    }
    if (ub < numeric_limits<uint16_t>::max()) {
        LOG("+ using 16bit DP table" << endl);
        return new KnapsackUpperBoundT<uint16_t>(inst, scale, max_bits);
    }
    if (ub < numeric_limits<uint32_t>::max()) {
        LOG("+ using 32bit DP table" << endl);
        return new KnapsackUpperBoundT<uint32_t>(inst, scale, max_bits);
    }
    if (ub < numeric_limits<uint64_t>::max()) {
        LOG("+ using 64bit DP table" << endl);
        return new KnapsackUpperBoundT<uint64_t>(inst, scale, max_bits);
    }
    return nullptr;
}
