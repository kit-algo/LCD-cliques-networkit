#include <algorithm>
#include <limits>
#include <unordered_map>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/unused.hpp>

#include <networkit/scd/TCE.hpp>

#include "LocalDegreeDirectedGraph.hpp"

namespace NetworKit {

TCE::TCE(const Graph &G, bool refine, bool useJaccard)
    : SelectiveCommunityDetector(G), refine(refine), useJaccard(useJaccard) {}

namespace {

#ifndef NDEBUG
/*
 * @deprecated Only used for debugging.
 */
double weightedEdgeScore(const Graph &G, node, node v, edgeweight uvWeight, double uDegree,
                         const std::unordered_map<node, double> &uNeighbors, bool useJaccard) {
    double nom = uvWeight;
    double vDegree = 0.0;

    G.forNeighborsOf(v, [&](node, node w, edgeweight vwWeight) {
        vDegree += vwWeight;
        auto uw = uNeighbors.find(w);
        if (uw != uNeighbors.end())
            nom += std::min(vwWeight, uw->second);
    });

    if (vDegree == 0.0)
        return 0.0;

    double denom = useJaccard ? (uDegree + vDegree - nom) : std::min(uDegree, vDegree);
    return nom / (denom * G.degree(v));
}
#endif

template <bool is_weighted>
std::set<node> expandSeedSet_internal(const Graph &G, const std::set<node> &s, bool refine,
                                      bool useJaccard) {
    // global data structures
    std::set<node> result = s;

    std::vector<double> weighted_degree;
    std::vector<double> node_score;

    std::vector<double> cutEdges;
    std::vector<bool> in_result;

    // heap that contains the nodes of the shell that still need to be considered
    class Compare {
        const std::vector<double> &similarity;

    public:
        Compare(const std::vector<double> &similarity) : similarity(similarity) {}

        bool operator()(node u, node v) { return similarity[u] > similarity[v]; }
    };

    tlx::d_ary_addressable_int_heap<node, 4, Compare> shell((Compare(node_score)));

    // data structures only for neighbors of a node
    std::vector<double> triangle_sum;

    auto node_added = [&](node, node, double weightedDegree) {
        weighted_degree.push_back(weightedDegree);
        node_score.push_back(0);
        cutEdges.push_back(0);

        in_result.push_back(false);
        triangle_sum.push_back(0);
    };

    LocalDegreeDirectedGraph<is_weighted, decltype(node_added)> localGraph(G, node_added);

    auto updateShell = [&](node u, node lu, bool noverify) -> double {
#ifdef NDEBUG
        tlx::unused(noverify); // only used in debug mode
#endif
        if (G.degree(u) == 0)
            return 0.0;

        double xDegree = weighted_degree[lu];

        localGraph.forTrianglesOf(
            u, [&](node lv, node y, double weight_uv, double weight_uy, double weight_vy) {
                if (is_weighted) {
                    triangle_sum[y] += std::min(weight_uv, weight_vy);
                    triangle_sum[lv] += std::min(weight_uy, weight_vy);
                } else {
                    triangle_sum[y] += 1;
                    triangle_sum[lv] += 1;
                }
            });

#ifndef NDEBUG
        std::unordered_map<node, double> uNeighbors;
        G.forNeighborsOf(
            u, [&](node, node v, edgeweight ew) { uNeighbors.insert(std::make_pair(v, ew)); });
#endif

        // collect counts and compute scores
        localGraph.forLocalNeighbors([&](node v, edgeweight weight) {
            if (!in_result[v]) {
                double nom = weight + triangle_sum[v];
                double wd = weighted_degree[v];

                double score_uv = 0.0;
                if (wd > 0.0) {
                    double denom = useJaccard ? (wd + xDegree - nom) : std::min(wd, xDegree);
                    score_uv = nom / (denom * G.degree(localGraph.toGlobal(v)));
                }

#ifndef NDEBUG
                assert(score_uv
                       == weightedEdgeScore(G, u, localGraph.toGlobal(v), weight, xDegree,
                                            uNeighbors, useJaccard));
#endif
                node_score[v] += score_uv;
                if (shell.contains(v)) {
                    shell.update(v);
                } else {
                    shell.push(v);
                }

                cutEdges[v] += weight;

#ifndef NDEBUG
                if (!noverify) {
                    double debugWeight = 0;
                    G.forNeighborsOf(localGraph.toGlobal(v), [&](node, node y, edgeweight ew) {
                        if (result.count(y)) {
                            debugWeight += ew;
                        }
                    });
                    assert(debugWeight == cutEdges[v]);
                }
#endif
            }

            // reset triangle_sum
            triangle_sum[v] = 0;
        });

        assert(
            std::all_of(triangle_sum.begin(), triangle_sum.end(), [](double s) { return s == 0; }));

        return xDegree;
    };

    // init community with seed set
    double volume = 0.0;
    double numCutEdges = 0.0;

    for (auto u : result) {
        node lu = localGraph.ensureNodeExists(u);
        in_result[lu] = true;
    }

    for (auto u : result) {
        volume += updateShell(u, localGraph.ensureNodeExists(u), true);
    }

    for (const auto &vv : cutEdges) {
        numCutEdges += vv;
    }

    // expand (main loop)
    while (!shell.empty()) {
        node uMax = shell.extract_top();
        node gUMax = localGraph.toGlobal(uMax);
        double uMaxVol = weighted_degree[uMax];

        double numCutEdgesNew = numCutEdges + uMaxVol - (2 * cutEdges[uMax]);
        double volumeNew = volume + uMaxVol;

#ifndef NDEBUG
        {
            double debugVolumeOld = 0.0;
            double debugCutsizeOld = 0.0;

            double debugVolumeNew = 0.0;
            double debugCutsizeNew = 0.0;

            for (auto u : result) {
                G.forEdgesOf(u, [&](node, node v, edgeweight ew) {
                    debugVolumeOld += ew;
                    debugVolumeNew += ew;
                    if (!result.count(v)) {
                        debugCutsizeOld += ew;

                        if (v != gUMax) {
                            debugCutsizeNew += ew;
                        }
                    }
                });
            }
            G.forEdgesOf(gUMax, [&](node, node v, edgeweight ew) {
                debugVolumeNew += ew;
                if (!result.count(v)) {
                    debugCutsizeNew += ew;
                }
            });

            assert(debugVolumeOld == volume);
            assert(debugVolumeNew == volumeNew);
            assert(debugCutsizeNew == numCutEdgesNew);
            assert(debugCutsizeOld == numCutEdges);
        }
#endif

        // if conductance decreases (i.e., improves)
        if ((numCutEdgesNew / volumeNew) < (numCutEdges / volume)) {
            result.insert(gUMax);
            in_result[uMax] = true;
            updateShell(gUMax, uMax, false);

            numCutEdges = numCutEdgesNew;
            volume = volumeNew;
        }
    }

    if (refine) {
        for (auto it = result.begin(); it != result.end();) {
            node u = *it;
            double uVol = 0;
            double uCutChange = 0;

            G.forNeighborsOf(u, [&](node, node v, edgeweight ew) {
                uVol += ew;
                if (result.count(v) > 0) {
                    uCutChange += ew;
                } else {
                    uCutChange -= ew;
                }
            });

            double numCutEdgesNew = numCutEdges + uCutChange;
            double volumeNew = volume - uVol;

#ifndef NDEBUG
            {
                double debugVolumeOld = 0.0;
                double debugCutsizeOld = 0.0;

                double debugVolumeNew = 0.0;
                double debugCutsizeNew = 0.0;

                for (auto v : result) {
                    G.forEdgesOf(v, [&](node, node x, edgeweight ew) {
                        debugVolumeOld += ew;
                        if (v != u) {
                            debugVolumeNew += ew;
                        }
                        if (!result.count(x)) {
                            debugCutsizeOld += ew;
                            if (v != u) {
                                debugCutsizeNew += ew;
                            }
                        }
                        if (x == u) {
                            debugCutsizeNew += ew;
                        }
                    });
                }

                assert(debugVolumeOld == volume);
                assert(debugVolumeNew == volumeNew);
                assert(debugCutsizeNew == numCutEdgesNew);
                assert(debugCutsizeOld == numCutEdges);
            }
#endif

            if ((numCutEdgesNew / volumeNew) < (numCutEdges / volume)) {
                numCutEdges = numCutEdgesNew;
                volume = volumeNew;

                TRACE("Removing ", u, " again");

                it = result.erase(it);
            } else {
                ++it;
            }
        }
    }

    return result;
}

} // namespace

std::set<node> TCE::expandOneCommunity(const std::set<node> &s) {
    if (G->isWeighted()) {
        return expandSeedSet_internal<true>(*G, s, refine, useJaccard);
    } else {
        return expandSeedSet_internal<false>(*G, s, refine, useJaccard);
    }
}

} /* namespace NetworKit */
