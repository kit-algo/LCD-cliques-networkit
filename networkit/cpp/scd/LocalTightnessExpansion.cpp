// networkit-format
#include <algorithm>
#include <limits>
#include <unordered_map>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/unused.hpp>

#include <networkit/scd/LocalTightnessExpansion.hpp>

#include "LocalDegreeDirectedGraph.hpp"

namespace NetworKit {

LocalTightnessExpansion::LocalTightnessExpansion(const Graph &g, double alpha)
    : SelectiveCommunityDetector(g), alpha(alpha) {}

namespace {

#ifndef NDEBUG
/*
 * @deprecated Only used for debugging.
 */
double weightedEdgeScore(const Graph &g, node u, node v, edgeweight ew, double uDegree,
                         const std::unordered_map<node, double> &uNeighbors) {
    double vDegree = 1.0; // w(u, u) = 1 per their definition
    double nom = ew;      // w(v, v) * w(v, u)

    g.forNeighborsOf(v, [&](node, node w, edgeweight vwWeight) {
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
    NodeAddedCallbackType nodeAddedCallback;
    std::vector<double> &triangleSum;

    void operator()(node u, node localId, double weightedDegree) {
        nodeAddedCallback(u, localId, weightedDegree);
        triangleSum.push_back(0);
    }
};

template <bool is_weighted, typename NodeAddedCallbackType>
class LocalGraph
    : public LocalDegreeDirectedGraph<is_weighted, InnerNodeAddedCallback<NodeAddedCallbackType>> {
private:
    // data structures only for neighbors of a node
    std::vector<double> triangleSum;

public:
    LocalGraph(const NetworKit::Graph &g, NodeAddedCallbackType nodeAddedCallback)
        : LocalDegreeDirectedGraph<is_weighted, InnerNodeAddedCallback<NodeAddedCallbackType>>(
            g, InnerNodeAddedCallback<NodeAddedCallbackType>{nodeAddedCallback, triangleSum}) {
        if (g.isWeighted() != is_weighted) {
            throw std::runtime_error("Error, weighted/unweighted status of input graph does not "
                                     "match is_weighted template parameter");
        }
    }

    template <typename F>
    void forNeighborsWithTrianglesOf(node u, F callback) {
        if (this->g.degree(u) == 0)
            return;

        this->forTrianglesOf(
            u, [&](node ln, edgeweight w) { triangleSum[ln] = 2 * w; },
            [&](node lv, node y, edgeweight weightUv, edgeweight weightUy, edgeweight weightVy) {
                if (is_weighted) {
                    triangleSum[y] += weightUv * weightVy;
                    triangleSum[lv] += weightUy * weightVy;
                } else {
                    triangleSum[y] += 1;
                    triangleSum[lv] += 1;
                }
            });

        this->forLocalNeighbors([&](node v, edgeweight w) { callback(v, w, triangleSum[v]); });
    }
};

template <bool is_weighted>
std::set<node> expandSeedSetInternal(const Graph &g, const std::set<node> &s, double alpha) {
    // global data structures
    // result community
    std::set<node> result;

    // This algorithm creates a local graph. It contains the community and its shell.
    // stores sqrt{sum_{v \in N(u)} w(u, v)^2} for every local node u
    std::vector<double> weightedDegree;

    // stores the sum of the similarities of all nodes to nodes in the community
    std::vector<double> nodeInternalSimilarity;
    std::vector<double> nodeExternalSimilarity;

    // indicates for every local node if it is in the community (true) or in the shell (false)
    std::vector<bool> inResult;

    std::vector<bool> inShell;

    // heap that contains the nodes of the shell that still need to be considered
    class Compare {
        const std::vector<double> &similarity;

    public:
        Compare(const std::vector<double> &similarity) : similarity(similarity) {}

        bool operator()(node u, node v) { return similarity[u] > similarity[v]; }
    };

    tlx::d_ary_addressable_int_heap<node, 4, Compare> shell((Compare(nodeInternalSimilarity)));

    double internalSimilarity = 0;
    double externalSimilarity = 0;

    auto addNode = [&](node u, node, double) {
        double wd = 1;
        if (is_weighted) {
            g.forNeighborsOf(u, [&](node, node, edgeweight w) { wd += w * w; });
        } else {
            wd += g.degree(u);
        }

        weightedDegree.push_back(std::sqrt(wd));

        nodeInternalSimilarity.push_back(.0);
        nodeExternalSimilarity.push_back(.0);
        inResult.push_back(false);
        inShell.push_back(false);
    };

    LocalGraph<is_weighted, decltype(addNode)> localGraph(g, addNode);

    auto updateShell = [&](node u, node lu) {
#ifndef NDEBUG
        std::unordered_map<node, double> uNeighbors;
        g.forNeighborsOf(
            u, [&](node, node v, edgeweight ew) { uNeighbors.insert(std::make_pair(v, ew)); });

        if (inShell[lu]) {
            double debugUExternalSimilarity = .0;
            double debugUInternalSimilarity = .0;
            localGraph.forNeighborsWithTrianglesOf(u, [&](node lv, edgeweight, double triangleSum) {
                double scoreUv = 0.0;
                double denom = weightedDegree[lv] * weightedDegree[lu];
                if (denom > 0.0) {
                    scoreUv = triangleSum / denom;
                }

                if (inResult[lv]) {
                    debugUInternalSimilarity += scoreUv;
                } else {
                    debugUExternalSimilarity += scoreUv;
                }
            });

            assert(std::abs(debugUExternalSimilarity - nodeExternalSimilarity[lu]) < 0.000001);
            assert(std::abs(debugUInternalSimilarity - nodeInternalSimilarity[lu]) < 0.000001);
        }
#endif

        std::vector<node> newShellNodes;

        localGraph.forNeighborsWithTrianglesOf(
            u, [&](node lv, edgeweight weight, double triangleSum) {
                // collect counts and compute scores
                double denom = weightedDegree[lv] * weightedDegree[lu];
                double scoreUv = triangleSum / denom;

#ifndef NDEBUG
                {
                    double debugScore = weightedEdgeScore(g, u, localGraph.toGlobal(lv), weight,
                                                          weightedDegree[lu], uNeighbors);
                    assert(scoreUv == debugScore);
                }
#else
            tlx::unused(weight); // needed for assert above
#endif

                nodeInternalSimilarity[lv] += scoreUv;

                if (inResult[lv]) {
                    externalSimilarity -= scoreUv;
                    internalSimilarity += 2 * scoreUv;
                    if (!inShell[lu]) {
                        nodeInternalSimilarity[lu] += scoreUv;
                    }
                    nodeExternalSimilarity[lv] -= scoreUv;
                } else {
                    externalSimilarity += scoreUv;

                    if (!inShell[lu]) {
                        nodeExternalSimilarity[lu] += scoreUv;
                    }

                    if (shell.contains(lv)) {
                        shell.update(lv);
                    } else {
                        shell.push(lv);
                    }

                    if (!inShell[lv]) {
                        inShell[lv] = true;
                        newShellNodes.push_back(lv);
                    } else {
                        nodeExternalSimilarity[lv] -= scoreUv;
                    }
                }
            });

        for (node s : newShellNodes) {
            localGraph.forNeighborsWithTrianglesOf(
                localGraph.toGlobal(s), [&](node lv, edgeweight, double triangleSum) {
                    if (!inResult[lv]) {
                        double scoreUv = 0.0;
                        double denom = weightedDegree[lv] * weightedDegree[s];
                        if (denom > 0.0) {
                            scoreUv = triangleSum / denom;
                        }

                        nodeExternalSimilarity[s] += scoreUv;
                    }
                });
        }
    };

    // init community with seed set
    for (node u : s) {
        node lu = localGraph.ensureNodeExists(u);
        result.insert(u);
        inResult[lu] = true;
        updateShell(u, lu);
    }

    // expand (main loop)
    while (!shell.empty()) {
        node uMax = shell.extract_top();
        node gUMax = localGraph.toGlobal(uMax);
#ifndef NDEBUG
        double internalSimilarityNew = internalSimilarity + 2 * nodeInternalSimilarity[uMax];
        double externalSimilarityNew =
            externalSimilarity - nodeInternalSimilarity[uMax] + nodeExternalSimilarity[uMax];

        {
            double debugInternalSimilarity = .0;
            double debugExternalSimilarity = .0;
            for (node u : result) {
                node lu = localGraph.ensureNodeExists(u);
                debugInternalSimilarity += nodeInternalSimilarity[lu];
                debugExternalSimilarity += nodeExternalSimilarity[lu];
            }

            assert(std::abs(debugInternalSimilarity - internalSimilarity) < 0.000001);
            assert(std::abs(debugExternalSimilarity - externalSimilarity) < 0.000001);
        }
#endif

        // if conductance decreases (i.e., improves)
        if (externalSimilarity / internalSimilarity
                - (alpha * nodeExternalSimilarity[uMax] - nodeInternalSimilarity[uMax])
                      / (2 * nodeInternalSimilarity[uMax])
            > 0) {
            result.insert(gUMax);
            inResult[uMax] = true;
            updateShell(gUMax, uMax);

            assert(std::abs(internalSimilarity - internalSimilarityNew) < 0.000001);
            assert(std::abs(externalSimilarity - externalSimilarityNew) < 0.000001);
        }
    }
    return result;
}

} // namespace

std::set<node> LocalTightnessExpansion::expandOneCommunity(const std::set<node> &s) {
    if (g->isWeighted()) {
        return expandSeedSetInternal<true>(*g, s, alpha);
    } else {
        return expandSeedSetInternal<false>(*g, s, alpha);
    }
}

} /* namespace NetworKit */
