#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#include <networkit/scd/LocalT.hpp>

#include "LocalDegreeDirectedGraph.hpp"

namespace NetworKit {

LocalT::LocalT(const Graph &G) : SelectiveCommunityDetector(G) {}

std::set<node> LocalT::expandOneCommunity(const std::set<node> &s) {
    // global data structures
    // result community
    std::set<node> result;

    // This algorithm creates a local graph. It contains the community and its shell.
    // stores the sum of the similarities of all nodes to nodes in the community
    std::vector<double> node_internal_triangles;
    std::vector<double> node_external_triangles;
    std::vector<double> node_semi_internal_triangles;

    // indicates for every local node if it is in the community (true) or in the shell (false)
    std::vector<bool> in_result;

    std::vector<bool> in_shell;

    count internal_triangles = 0, external_triangles = 0;

    auto addNode = [&](node, node, double) {
        node_internal_triangles.push_back(0);
        node_external_triangles.push_back(0);
        node_semi_internal_triangles.push_back(0);
        in_result.push_back(false);
        in_shell.push_back(false);
    };

    LocalDegreeDirectedGraph<false, decltype(addNode)> local_graph(*G, addNode);

    std::unordered_set<node> shell;

    auto updateShell = [&](node u) {
        local_graph.forTrianglesOf(u, [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
            // collect counts and compute scores
            // new completely internal triangle, was not counted previously
            if (in_result[lv] && in_result[lw]) {
                ++node_internal_triangles[lv];
                ++node_internal_triangles[lw];
                ++internal_triangles;
            } else if (in_result[lv] || in_result[lw]) {
                --external_triangles;

                if (in_result[lv]) {
                    ++node_internal_triangles[lw];
                    --node_semi_internal_triangles[lw];
                } else {
                    ++node_internal_triangles[lv];
                    --node_semi_internal_triangles[lv];
                }
            } else {
                ++external_triangles;
                if (in_shell[lv]) {
                    --node_external_triangles[lv];
                }
                ++node_semi_internal_triangles[lv];

                if (in_shell[lw]) {
                    --node_external_triangles[lw];
                }
                ++node_semi_internal_triangles[lw];
            }
        });

        G->forNeighborsOf(u, [&](node v) {
            node s = local_graph.ensureNodeExists(v);

            if (!in_shell[s] && !in_result[s]) {
                shell.insert(s);
                in_shell[s] = true;

                local_graph.forTrianglesOf(
                    local_graph.toGlobal(s),
                    [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
                        if (!in_result[lv] && !in_result[lw]) {
                            ++node_external_triangles[s];
                        }
                    });
            }
        });

#ifndef NDEBUG
        count debug_external_triangles = 0;
        count debug_internal_triangles = 0;
        for (node u : result) {
            count u_external = 0, u_internal = 0;
            local_graph.forTrianglesOf(u,
                                       [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
                                           if (in_result[lv] && in_result[lw]) {
                                               ++u_internal;
                                           } else if (!in_result[lv] && !in_result[lw]) {
                                               ++u_external;
                                           }
                                       });

            node lu = local_graph.ensureNodeExists(u);

            assert(u_internal == node_internal_triangles[lu]);
            assert(u_internal == node_internal_triangles[lu]);

            debug_external_triangles += u_external;
            debug_internal_triangles += u_internal;
        }

        assert(debug_external_triangles == external_triangles);
        assert(debug_internal_triangles / 3 == internal_triangles);

        for (node ls : shell) {
            count s_external = 0, s_internal = 0, s_semi_internal = 0;
            local_graph.forTrianglesOf(local_graph.toGlobal(ls),
                                       [&](node lv, node lw, edgeweight, edgeweight, edgeweight) {
                                           if (in_result[lv] && in_result[lw]) {
                                               ++s_internal;
                                           } else if (in_result[lv] || in_result[lw]) {
                                               ++s_semi_internal;
                                           } else {
                                               ++s_external;
                                           }
                                       });

            assert(s_semi_internal == node_semi_internal_triangles[ls]);
            assert(s_internal == node_internal_triangles[ls]);
            assert(s_external == node_external_triangles[ls]);
        }
#endif
    };

    // init community with seed set
    for (node u : s) {
        node lu = local_graph.ensureNodeExists(u);
        result.insert(u);
        in_result[lu] = true;
        updateShell(u);
    }

    auto get_score = [](count int_triangles, count ext_triangles) -> count {
        return std::max<int64_t>(
            0, int_triangles
                   * (static_cast<int64_t>(int_triangles) - static_cast<int64_t>(ext_triangles)));
    };

    // expand (main loop)
    node uMax = none;
    do {

        uMax = none;
        count best_score = get_score(internal_triangles, external_triangles);
        count best_external_triangles = none;

        for (node lv : shell) {
            count new_internal_triangles = internal_triangles + node_internal_triangles[lv];
            count new_external_triangles =
                external_triangles + node_external_triangles[lv] - node_semi_internal_triangles[lv];

            count new_score = get_score(new_internal_triangles, new_external_triangles);

            if (new_score > best_score
                || (new_score == best_score && new_external_triangles < best_external_triangles)) {
                uMax = lv;
                best_score = new_score;
                best_external_triangles = new_external_triangles;
            }
        }

        if (uMax != none) {
            node gUMax = local_graph.toGlobal(uMax);
            result.insert(gUMax);
            shell.erase(uMax);
            in_result[uMax] = true;
            updateShell(gUMax);
            assert(external_triangles == best_external_triangles);
            assert(get_score(internal_triangles, external_triangles) == best_score);
        }

    } while (uMax != none);
    return result;
}

} /* namespace NetworKit */
