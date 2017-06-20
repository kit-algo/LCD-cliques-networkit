#include <algorithm>
#include <limits>
#include <unordered_map>

#include <tlx/container/d_ary_addressable_int_heap.hpp>
#include <tlx/unused.hpp>

#include <networkit/scd/TCE.hpp>

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

    std::vector<node> local_to_global_id;
    std::unordered_map<node, node> global_to_local_id;

    std::vector<node> head;
    std::vector<double> head_weight;

    struct local_node {
        std::size_t first_head;
        std::size_t last_head;
    };

    std::vector<local_node> head_info;
    std::vector<count> degree;

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
    std::vector<bool> is_neighbor;
    std::vector<double> triangle_sum;
    std::vector<double> neighbor_weight;

    auto ensureNodeExists = [&](node u) -> node {
        auto gtl_it = global_to_local_id.find(u);

        if (gtl_it == global_to_local_id.end()) {
            node local_id = local_to_global_id.size();
            local_to_global_id.push_back(u);
            global_to_local_id[u] = local_id;

            double weightedDegree = 0;

            size_t my_degree = 0;

            index my_begin = head.size();
            index my_end = my_begin;
            index nh = my_end; // next head if all potential out neighbors were inserted

            G.forEdgesOf(u, [&](node, node v, edgeweight weight) {
                auto it = global_to_local_id.find(v);

                if (it != global_to_local_id.end()) {
                    ++my_degree;

                    auto &n_degree = degree[it->second];
                    ++n_degree;

                    head.push_back(it->second);
                    if (is_weighted) {
                        head_weight.push_back(weight);
                    }

                    ++my_end;
                }

                weightedDegree += weight;
                ++nh;
            });

            assert(my_end == head.size());
            assert(!is_weighted || my_end == head_weight.size());
            assert(nh >= my_end);

            head.resize(nh);
            if (is_weighted) {
                head_weight.resize(nh);
            }

            degree.push_back(my_degree);

            for (index i = my_begin; i < my_end;) {
                node ln = head[i];
                edgeweight weight = is_weighted ? head_weight[i] : 1;

                auto &n_degree = degree[ln];
                auto &n_info = head_info[ln];

                // before pushing into n_edges check if we can get rid of any neighbors.
                for (index ni = n_info.first_head; ni < n_info.last_head;) {
                    node nn = head[ni];
                    if (n_degree <= degree[nn]) {
                        ++ni;
                    } else {
                        head[ni] = head[--n_info.last_head];

                        if (is_weighted) {
                            head_weight[head_info[nn].last_head] = head_weight[ni];
                            head_weight[ni] = head_weight[n_info.last_head];
                        }

                        head[head_info[nn].last_head++] = ln;
                    }
                }

                if (my_degree < n_degree) {
                    ++i;
                } else {
                    if (is_weighted) {
                        head_weight[n_info.last_head] = weight;
                    }

                    head[n_info.last_head++] = local_id;

                    --my_end;
                    head[i] = head[my_end];
                    if (is_weighted) {
                        head_weight[i] = head_weight[my_end];
                    }

                    // assert(head_info[ln].last_head <= head_info[ln+1].first_head);
                }
            }

            head_info.emplace_back(local_node{my_begin, my_end});

            // head_info.back().last_head = head.size();

            weighted_degree.push_back(weightedDegree);
            node_score.push_back(0);
            cutEdges.push_back(0);

            is_neighbor.push_back(false);
            in_result.push_back(false);
            triangle_sum.push_back(0);

            if (is_weighted) {
                neighbor_weight.push_back(0);
            }

            return local_id;
        } else {
            return gtl_it->second;
        }
    };

    std::vector<node> current_local_neighbor_ids;

    auto updateShell = [&](node u, bool noverify) -> double {
#ifdef NDEBUG
        tlx::unused(noverify); // only used in debug mode
#endif

        if (G.degree(u) == 0.0)
            return 0.0;

        double xDegree = 0.0;

        current_local_neighbor_ids.clear();
        current_local_neighbor_ids.reserve(G.degree(u));

        G.forNeighborsOf(u, [&](node, node v, edgeweight weight) {
            xDegree += weight;

            node local_id = ensureNodeExists(v);
            current_local_neighbor_ids.emplace_back(local_id);
            is_neighbor[local_id] = true;
            if (is_weighted) {
                neighbor_weight[local_id] = weight;
            }
            triangle_sum[local_id] = 0;
        });

        for (node lv : current_local_neighbor_ids) {
            // count (weighted) triangles using only out-neighbors
            double triangles_lv = 0;
            for (index ti = head_info[lv].first_head; ti < head_info[lv].last_head; ++ti) {
                node y = head[ti];

                if (is_neighbor[y]) {
                    if (is_weighted) {
                        triangle_sum[y] += std::min(neighbor_weight[lv], head_weight[ti]);
                        triangles_lv += std::min(neighbor_weight[y], head_weight[ti]);
                    } else {
                        triangle_sum[y] += 1;
                        triangles_lv += 1;
                    }
                }
            }

            triangle_sum[lv] += triangles_lv;
        }
#ifndef NDEBUG
        std::unordered_map<node, double> uNeighbors;
        G.forNeighborsOf(
            u, [&](node, node v, edgeweight ew) { uNeighbors.insert(std::make_pair(v, ew)); });
#endif

        // collect counts and compute scores
        for (node v : current_local_neighbor_ids) {
            if (!in_result[v]) {

                edgeweight weight = (is_weighted ? neighbor_weight[v] : 1);
                double nom = weight + triangle_sum[v];
                double wd = weighted_degree[v];

                double score_uv = 0.0;
                if (wd > 0.0) {
                    double denom = useJaccard ? (wd + xDegree - nom) : std::min(wd, xDegree);
                    score_uv = nom / (denom * G.degree(local_to_global_id[v]));
                }

#ifndef NDEBUG
                assert(score_uv
                       == weightedEdgeScore(G, u, local_to_global_id[v], weight, xDegree,
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
                    G.forNeighborsOf(local_to_global_id[v], [&](node, node y, edgeweight ew) {
                        if (result.count(y)) {
                            debugWeight += ew;
                        }
                    });
                    assert(debugWeight == cutEdges[v]);
                }
#endif
            }
        }

        for (node v : current_local_neighbor_ids) {
            is_neighbor[v] = false;
        }

        return xDegree;
    };

    // init community with seed set
    double volume = 0.0;
    double numCutEdges = 0.0;

    for (auto u : result) {
        node lu = ensureNodeExists(u);
        in_result[lu] = true;
    }

    for (auto u : result) {
        volume += updateShell(u, true);
    }

    for (const auto &vv : cutEdges) {
        numCutEdges += vv;
    }

    // expand (main loop)
    while (!shell.empty()) {
        node uMax = shell.extract_top();
        node gUMax = local_to_global_id[uMax];
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
            updateShell(gUMax, false);

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
