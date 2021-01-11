// networkit-format
#include <algorithm>
#include <limits>
#include <unordered_map>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/unused.hpp>

#include <networkit/scd/LocalTightnessExpansion.hpp>

#include "LocalDegreeDirectedGraph.hpp"

namespace NetworKit {

LocalTightnessExpansion::LocalTightnessExpansion(const Graph &G, double alpha)
    : SelectiveCommunityDetector(G), alpha(alpha) {}

namespace {

#ifndef NDEBUG
/*
 * @deprecated Only used for debugging.
 */
double weightedEdgeScore(const Graph &G, node u, node v, edgeweight ew, double uDegree,
                         const std::unordered_map<node, double> &uNeighbors) {
    double vDegree = 1.0; // w(u, u) = 1 per their definition
    double nom = ew;      // w(v, v) * w(v, u)

    G.forNeighborsOf(v, [&](node, node w, edgeweight vwWeight) {
        vDegree += vwWeight * vwWeight;
        auto uw = uNeighbors.find(w);
        if (uw != uNeighbors.end())
            nom += vwWeight * uw->second;
        else if (w == u)
            nom += vwWeight; // w(u, u) * w(u, v)
    });

    vDegree = std::sqrt(vDegree);

    double denom = vDegree * uDegree;

    return nom / denom;
}
#endif

template <typename NodeAddedCallbackType>
struct InnerNodeAddedCallback {
    NodeAddedCallbackType node_added_callback;
    std::vector<double> &triangle_sum;

    void operator()(node u, node local_id, double weightedDegree) {
        node_added_callback(u, local_id, weightedDegree);
        triangle_sum.push_back(0);
    }
};

template <bool is_weighted, typename NodeAddedCallbackType>
class LocalGraph
    : public LocalDegreeDirectedGraph<is_weighted, InnerNodeAddedCallback<NodeAddedCallbackType>> {
private:
    // data structures only for neighbors of a node
    std::vector<double> triangle_sum;

public:
    LocalGraph(const NetworKit::Graph &G, NodeAddedCallbackType node_added_callback)
        : LocalDegreeDirectedGraph<is_weighted, InnerNodeAddedCallback<NodeAddedCallbackType>>(
            G, InnerNodeAddedCallback<NodeAddedCallbackType>{node_added_callback, triangle_sum}) {
        if (G.isWeighted() != is_weighted) {
            throw std::runtime_error("Error, weighted/unweighted status of input graph does not "
                                     "match is_weighted template parameter");
        }
    }

    template <typename F>
    void forNeighborsWithTrianglesOf(node u, F callback) {
        if (this->G.degree(u) == 0)
            return;

        this->forTrianglesOf(
            u, [&](node ln, edgeweight w) { triangle_sum[ln] = 2 * w; },
            [&](node lv, node y, edgeweight weight_uv, edgeweight weight_uy, edgeweight weight_vy) {
                if (is_weighted) {
                    triangle_sum[y] += weight_uv * weight_vy;
                    triangle_sum[lv] += weight_uy * weight_vy;
                } else {
                    triangle_sum[y] += 1;
                    triangle_sum[lv] += 1;
                }
            });

        this->forLocalNeighbors([&](node v, edgeweight w) { callback(v, w, triangle_sum[v]); });
    }
};

template <bool is_weighted>
std::set<node> expandSeedSet_internal(const Graph &G, const std::set<node> &s, double alpha) {
    // global data structures
    // result community
    std::set<node> result;

    // This algorithm creates a local graph. It contains the community and its shell.
    // stores sqrt{sum_{v \in N(u)} w(u, v)^2} for every local node u
    std::vector<double> weighted_degree;

    // stores the sum of the similarities of all nodes to nodes in the community
    std::vector<double> node_internal_similarity;
    std::vector<double> node_external_similarity;

    // indicates for every local node if it is in the community (true) or in the shell (false)
    std::vector<bool> in_result;

    std::vector<bool> in_shell;

    // heap that contains the nodes of the shell that still need to be considered
    class Compare {
        const std::vector<double> &similarity;

    public:
        Compare(const std::vector<double> &similarity) : similarity(similarity) {}

        bool operator()(node u, node v) { return similarity[u] > similarity[v]; }
    };

    tlx::d_ary_addressable_int_heap<node, 4, Compare> shell((Compare(node_internal_similarity)));

    double internal_similarity = 0;
    double external_similarity = 0;

    auto addNode = [&](node u, node, double) {
        double wd = 1;
        if (is_weighted) {
            G.forNeighborsOf(u, [&](node, node, edgeweight w) { wd += w * w; });
        } else {
            wd += G.degree(u);
        }

        weighted_degree.push_back(std::sqrt(wd));

        node_internal_similarity.push_back(.0);
        node_external_similarity.push_back(.0);
        in_result.push_back(false);
        in_shell.push_back(false);
    };

    LocalGraph<is_weighted, decltype(addNode)> local_graph(G, addNode);

    auto updateShell = [&](node u, node lu) {
#ifndef NDEBUG
        std::unordered_map<node, double> uNeighbors;
        G.forNeighborsOf(
            u, [&](node, node v, edgeweight ew) { uNeighbors.insert(std::make_pair(v, ew)); });

        if (in_shell[lu]) {
            double debug_u_external_similarity = .0;
            double debug_u_internal_similarity = .0;
            local_graph.forNeighborsWithTrianglesOf(
                u, [&](node lv, edgeweight, double triangle_sum) {
                    double score_uv = 0.0;
                    double denom = weighted_degree[lv] * weighted_degree[lu];
                    if (denom > 0.0) {
                        score_uv = triangle_sum / denom;
                    }

                    if (in_result[lv]) {
                        debug_u_internal_similarity += score_uv;
                    } else {
                        debug_u_external_similarity += score_uv;
                    }
                });

            assert(std::abs(debug_u_external_similarity - node_external_similarity[lu]) < 0.000001);
            assert(std::abs(debug_u_internal_similarity - node_internal_similarity[lu]) < 0.000001);
        }
#endif

        std::vector<node> new_shell_nodes;

        local_graph.forNeighborsWithTrianglesOf(
            u, [&](node lv, edgeweight weight, double triangle_sum) {
                // collect counts and compute scores
                double denom = weighted_degree[lv] * weighted_degree[lu];
                double score_uv = triangle_sum / denom;

#ifndef NDEBUG
                {
                    double debug_score = weightedEdgeScore(G, u, local_graph.toGlobal(lv), weight,
                                                           weighted_degree[lu], uNeighbors);
                    assert(score_uv == debug_score);
                }
#else
            tlx::unused(weight); // needed for assert above
#endif

                node_internal_similarity[lv] += score_uv;

                if (in_result[lv]) {
                    external_similarity -= score_uv;
                    internal_similarity += 2 * score_uv;
                    if (!in_shell[lu]) {
                        node_internal_similarity[lu] += score_uv;
                    }
                    node_external_similarity[lv] -= score_uv;
                } else {
                    external_similarity += score_uv;

                    if (!in_shell[lu]) {
                        node_external_similarity[lu] += score_uv;
                    }

                    if (shell.contains(lv)) {
                        shell.update(lv);
                    } else {
                        shell.push(lv);
                    }

                    if (!in_shell[lv]) {
                        in_shell[lv] = true;
                        new_shell_nodes.push_back(lv);
                    } else {
                        node_external_similarity[lv] -= score_uv;
                    }
                }
            });

        for (node s : new_shell_nodes) {
            local_graph.forNeighborsWithTrianglesOf(
                local_graph.toGlobal(s), [&](node lv, edgeweight, double triangle_sum) {
                    if (!in_result[lv]) {
                        double score_uv = 0.0;
                        double denom = weighted_degree[lv] * weighted_degree[s];
                        if (denom > 0.0) {
                            score_uv = triangle_sum / denom;
                        }

                        node_external_similarity[s] += score_uv;
                    }
                });
        }
    };

    // init community with seed set
    for (node u : s) {
        node lu = local_graph.ensureNodeExists(u);
        result.insert(u);
        in_result[lu] = true;
        updateShell(u, lu);
    }

    // expand (main loop)
    while (!shell.empty()) {
        node uMax = shell.extract_top();
        node gUMax = local_graph.toGlobal(uMax);
#ifndef NDEBUG
        double internal_similarity_new = internal_similarity + 2 * node_internal_similarity[uMax];
        double external_similarity_new =
            external_similarity - node_internal_similarity[uMax] + node_external_similarity[uMax];

        {
            double debug_internal_similarity = .0;
            double debug_external_similarity = .0;
            for (node u : result) {
                node lu = local_graph.ensureNodeExists(u);
                debug_internal_similarity += node_internal_similarity[lu];
                debug_external_similarity += node_external_similarity[lu];
            }

            assert(std::abs(debug_internal_similarity - internal_similarity) < 0.000001);
            assert(std::abs(debug_external_similarity - external_similarity) < 0.000001);
        }
#endif

        // if conductance decreases (i.e., improves)
        if (external_similarity / internal_similarity
                - (alpha * node_external_similarity[uMax] - node_internal_similarity[uMax])
                      / (2 * node_internal_similarity[uMax])
            > 0) {
            result.insert(gUMax);
            in_result[uMax] = true;
            updateShell(gUMax, uMax);

            assert(std::abs(internal_similarity - internal_similarity_new) < 0.000001);
            assert(std::abs(external_similarity - external_similarity_new) < 0.000001);
        }
    }
    return result;
}

} // namespace

std::set<node> LocalTightnessExpansion::expandOneCommunity(const std::set<node> &s) {
    if (G->isWeighted()) {
        return expandSeedSet_internal<true>(*G, s, alpha);
    } else {
        return expandSeedSet_internal<false>(*G, s, alpha);
    }
}

} /* namespace NetworKit */
