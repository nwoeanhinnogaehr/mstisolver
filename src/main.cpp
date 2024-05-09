#include <ctime>
#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;
using namespace std::chrono;

// awful global hack, so that we can get the node count even on timeout
uint64_t n_nodes = 0, n_pruned = 0;

#include "instance.hpp"
#include "knapsack_ub.hpp"
#include "bnb.hpp"
#include "bnb_blocker.hpp"
#include "greedy.hpp"
#include "output.hpp"

#include <clipp.h>
using namespace clipp;

void print_node_count(int sig) {
    OUT("num_bnb_nodes " << n_nodes << endl);
    OUT("num_bnb_pruned " << n_pruned << endl);
    signal(sig, SIG_DFL);
    kill(getpid(), sig);
}

int main(int argc, char **argv) {
    auto start_clock = std::clock();
    auto start_time = high_resolution_clock::now();

    bool quiet = false, use_file = false, write_canonical = false,
         write_instance = false,
         use_dp_ub = true, use_greedy = true, gen_inst = false,
         use_canonicalize = true, use_binsearch = false;
    int gen_n, gen_m, gen_max_lo_wt, gen_max_up_wt, rng_seed = -1;
    int cap_override = -1, r_override = -1;
    double cap_override_fraction = -1, r_override_fraction = -1;
    double gen_cap_scale;
    int64_t dp_cost_scale = 1;
    int dp_prefix_bits = 30;
    string filename, canonical_filename, write_inst_filename;

    auto cli =
        (((option("-f").set(use_file) & value("filename", filename)) %
              "specify instance to solve",
          (option("-cap") & value("cap", cap_override)) %
              "override instance capacity",
          (option("-cap-fraction") &
           value("cap-fraction", cap_override_fraction)) %
              "override instance capacity (specified by fraction of mincut - "
              "1)",
          (option("-r") & value("target_weight", r_override)) %
              "override instance r-value (target_weight)",
          (option("-r-fraction") &
           value("target_weight-fraction", r_override_fraction)) %
              "override instance r-value (specified by fraction)",
          (option("-g").set(gen_inst) %
               "generate a MST interdiction instance, with the following "
               "parameters:" &
           integer("n", gen_n) % "number of vertices" &
           integer("m", gen_m) % "number of edges" &
           integer("max-lo-wt", gen_max_lo_wt) % "maximum MST edge weight" &
           integer("max-up-wt", gen_max_up_wt) % "maximum interdiction cost" &
           value("cap-scale", gen_cap_scale) %
               "how much to scale down the interdiction budget (1.0 = maximum "
               "possible, 0.0 = no budget)") %
              "instance generation",
          (option("-write-inst").set(write_instance) &
           value("filename", write_inst_filename)) %
              "write original instance to a file",
          (option("-write-canon").set(write_canonical) &
           value("filename", canonical_filename)) %
              "write canonicalized instance to a file",
          (option("-q").set(quiet)) % "quiet mode: do not log to stderr") %
             "input / output",
         ((option("-no-canonicalize").set(use_canonicalize, false) %
           "disable instance canonicalization"),
          (option("-use-binsearch").set(use_binsearch, true) %
           "use binary search to solve minimum edge blocker instances"),
          ((option("-no-greedy").set(use_greedy, false)) %
               "disable greedy lower bound") %
              "greedy lower bound",
          ((option("-no-dp-ub").set(use_dp_ub, false)) %
               "disable DP upper bound",
           (option("-dp-cost-scale") & value("scale", dp_cost_scale)) %
               "scale down costs and capacity for DP bound",
           (option("-max-dp-prefix-bits") & value("bits", dp_prefix_bits)) %
               "maximum number of DP prefix bits to compute") %
              "dynamic programming upper bound") %
             "new algorithms",
         (option("-seed") & integer("seed", rng_seed)) %
             "seed to use for random number generation");

    if (!parse(argc, argv, cli)) {
        OUT(make_man_page(cli, argv[0]));
        return 0;
    }
    if (quiet) _log = false;

    // read/generate instance
    Instance inst;
    if (use_file) {
        ifstream inst_file(filename);
        ASSERT(inst_file.good());
        inst = Instance::read_from(inst_file);
        inst_file.close();
    } else if (gen_inst) {
        inst =
            Instance::generate_random(gen_n, gen_m, gen_max_lo_wt,
                                      gen_max_up_wt, gen_cap_scale, rng_seed);
    } else {
        LOG("+ no file or generator specified, reading instance from stdin" << endl);
        inst = Instance::read_from(cin);
    }

    // canonicalization can change the max spanning tree, so we do this first
    int64_t max_st = 0;
    UnionFind uf_max(inst.n());
    for (int i = inst.m() - 1; i >= 0; i--) {
        Edge &e = inst.edge(i);
        int a = uf_max.find(e.s);
        int b = uf_max.find(e.t);
        if (a != b) {
            max_st += e.lo_wt;
            uf_max.join(a, b);
        }
    }
    OUT("max_spanning_tree " << max_st << endl);
    OUT("min_cut " << inst.min_cut() << endl);
    Soln initial_mst(inst);

    // handle overrides
    if (r_override >= 0) {
        inst.set_r(r_override);
        LOG("+ set r to " << inst.r() << endl);
    }
    if (r_override_fraction >= 0) {
        inst.set_r(r_override_fraction * (max_st - initial_mst.weight) + initial_mst.weight);
        LOG("+ set r to " << inst.r() << endl);
    }
    if (cap_override >= 0)
        inst.set_cap(cap_override);
    if (cap_override_fraction >= 0)
        ASSERT(false); //inst.set_cap(cap);

    // check instance validity
    if (inst.problem_type() == ProblemType::MSTInterdiction) {
        if (inst.cap() < 0) {
            LOG("instance has no capacity, please specify -cap or -cap-fraction" << endl);
            OUT("status err" << endl);
            return 1;
        }
        if (!inst.is_valid()) {
            LOG("instance is invalid, leader can cut the graph." << endl);
            OUT("status err" << endl);
            return 1;
        }
    } else {
        if (inst.r() < 0) {
            LOG("instance has no r value (mst weight target), please specify -r or -r-fraction" << endl);
            OUT("status err" << endl);
            return 1;
        }
    }

    // write instance
    if (write_instance) {
        ofstream inst_out_file(write_inst_filename);
        ASSERT(inst_out_file.good());
        inst.write_to(inst_out_file);
        inst_out_file.close();
    }

    if (inst.min_cut() == 0) {
        LOG("+ instance is trivial: graph can be disconnected for free" << endl);
        OUT("sol_val cut" << endl);
        OUT("status ok" << endl);
        auto end_clock = std::clock();
        auto end_time = high_resolution_clock::now();
        OUT("wallclock_time " << duration<double, std::milli>(end_time - start_time).count() << std::endl);
        OUT("cpu_time " << 1000.0 * (end_clock - start_clock) / CLOCKS_PER_SEC << std::endl);
        print_node_count(SIGINT);
        return 0;
    }

    if (inst.problem_type() == ProblemType::MSTInterdiction) {
        // canonicalize
        Instance orig = inst;
        if (use_canonicalize)
            inst = inst.canonicalize();
        ASSERT(inst.is_valid()); // If this fails, canonicalization is definitely broken
        if (write_canonical) {
            ofstream inst_out_file(canonical_filename);
            ASSERT(inst_out_file.good());
            inst.write_to(inst_out_file);
            inst_out_file.close();
        }

        // compute MST
        Soln mst(inst);
        Soln initial_lb = mst;
        OUT("min_spanning_tree " << initial_lb.weight << endl);

        // compute dynamic programming upper bound
        KnapsackUpperBound *kp_ub = nullptr;
        if (use_dp_ub) {
            kp_ub = make_knapsack_ub(inst, dp_cost_scale, dp_prefix_bits+1);
            kp_ub->compute_in_background();
        }

        // compute greedy lower bound
        if (use_greedy) {
            Greedy greedy(inst);
            initial_lb = greedy.compute();
            OUT("greedy_lower_bound " << initial_lb.weight << endl);
        }

        Soln sol = mst;
        // solve with branch-and-bound
        BranchAndBound bnb(inst, initial_lb, kp_ub);
        signal(SIGTERM, print_node_count);
        signal(SIGINT, print_node_count);
        sol = bnb.solve();
        delete kp_ub; // kill background thread

        ASSERT(sol.weight >= initial_lb.weight);
        ASSERT(sol.weight <= max_st);
        ASSERT(mst.weight <= sol.weight);
        OUT("sol_val " << sol.weight << endl);
        OUT("X ";
        for (int i = 0; i < inst.m(); i++)
            cout << sol.x[i];
        cout << endl << "Y ";
        for (int i = 0; i < inst.m(); i++) {
            cout << sol.y[i];
        }
        cout << endl);
    } else { // minimum edge blocker instance
        if (!use_binsearch) {
            // canonicalize
            Instance orig = inst;
            if (use_canonicalize) {
                int64_t r = inst.r();
                inst = Instance(inst.n(), inst.m(), inst.min_cut()-1, vector<Edge>(inst.edges().begin(), inst.edges().end()), -1, inst.min_cut());
                inst = inst.canonicalize();
                ASSERT(inst.is_valid()); // If this fails, canonicalization is definitely broken
                inst.set_r(r);
            }
            if (write_canonical) {
                ofstream inst_out_file(canonical_filename);
                ASSERT(inst_out_file.good());
                inst.write_to(inst_out_file);
                inst_out_file.close();
            }

            // compute MST
            Soln mst(inst);
            Soln initial_lb = mst;
            OUT("min_spanning_tree " << initial_lb.weight << endl);

            // compute dynamic programming upper bound
            KnapsackUpperBound *kp_ub = nullptr;
            if (use_dp_ub) {
                kp_ub = make_knapsack_ub(inst, dp_cost_scale, dp_prefix_bits+1);
                kp_ub->compute_in_background();
            }

            Soln sol = mst;
            // solve with branch-and-bound
            BlockerBranchAndBound bnb(inst, kp_ub);
            signal(SIGTERM, print_node_count);
            signal(SIGINT, print_node_count);
            sol = bnb.solve();
            delete kp_ub; // kill background thread

            int64_t cost = 0;
            for (int i = 0; i < inst.m(); i++)
                if (sol.x[i]) cost += inst.edge(i).up_wt;
            if (cost == 0) {
                OUT("sol_val cut" << endl);
            } else {
                OUT("X ";
                for (int i = 0; i < inst.m(); i++)
                    cout << sol.x[i];
                cout << endl << "Y ";
                for (int i = 0; i < inst.m(); i++)
                    cout << sol.y[i];
                cout << endl;
                cout << "sol_val " << cost << endl);
            }
        } else {
            signal(SIGTERM, print_node_count);
            signal(SIGINT, print_node_count);
        
            // compute MST
            Soln mst(inst);
            OUT("min_spanning_tree " << mst.weight << endl);

            // compute dynamic programming upper bound
            KnapsackUpperBound *kp_ub = nullptr;
            if (use_dp_ub) {
                inst.set_cap(max(0l, inst.min_cut()-1));
                if (use_canonicalize)
                    kp_ub = make_knapsack_ub(inst.canonicalize(), dp_cost_scale, dp_prefix_bits+1);
                else
                    kp_ub = make_knapsack_ub(inst, dp_cost_scale, dp_prefix_bits+1);
                kp_ub->compute_in_background();
            }

            // binary search to find the optimal capacity
            int64_t lo = 0, hi = inst.min_cut();
            while (lo < hi) {
                int64_t mid = (lo + hi)/2;
                ASSERT(mid < inst.min_cut());
                LOG("+ binary search, current capacity " << mid << endl);

                if (kp_ub && kp_ub->prefix_bits_computed() >= 0) {
                    int64_t upper_bound = kp_ub->get_bound(0, mid, 0) + mst.weight;
                    LOG("+ upper bound: " << upper_bound << endl);
                    if (upper_bound < inst.r()) {
                        LOG("+ r exceeds upper bound, need more capacity" << endl);
                        lo = mid + 1;
                        continue;
                    }
                }

                Soln initial_lb = mst;

                // canonicalize
                Instance msti_inst(inst.n(), inst.m(), mid, vector<Edge>(inst.edges().begin(), inst.edges().end()), -1, inst.min_cut());
                if (use_canonicalize)
                    msti_inst = msti_inst.canonicalize();
                ASSERT(msti_inst.is_valid()); // If this fails, canonicalization is definitely broken

                // compute greedy lower bound
                if (use_greedy) {
                    Greedy greedy(msti_inst);
                    initial_lb = greedy.compute();
                }

                Soln sol = mst;
                // solve with branch-and-bound
                BranchAndBound bnb(msti_inst, initial_lb, kp_ub, inst.r()-1);
                sol = bnb.solve();
                LOG("+ nodes: " << n_nodes << endl);
                ASSERT(sol.weight >= initial_lb.weight);
                ASSERT(sol.weight <= max_st);
                ASSERT(mst.weight <= sol.weight);
                LOG("+ result: " << sol.weight << endl);
                if (sol.weight >= inst.r()) {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            if (lo == inst.min_cut()) {
                OUT("sol_val cut" << endl);
            } else {
                OUT("sol_val " << lo << endl);
            }
        }
    }
    OUT("status ok" << endl);
    auto end_clock = std::clock();
    auto end_time = high_resolution_clock::now();
    OUT("wallclock_time " << duration<double, std::milli>(end_time - start_time).count() << endl);
    OUT("cpu_time " << 1000.0 * (end_clock - start_clock) / CLOCKS_PER_SEC << endl);
    print_node_count(SIGINT);
}
