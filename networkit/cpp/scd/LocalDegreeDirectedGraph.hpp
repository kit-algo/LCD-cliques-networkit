#ifndef NETWORKIT_SCD_LOCALDEGREEDIRECTEDGRAPH_HPP
#define NETWORKIT_SCD_LOCALDEGREEDIRECTEDGRAPH_HPP

#include <unordered_map>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Graph for local community detection that stores edges in directed
 * form in just one direction depending on the node's degrees and
 * thus allows efficient triangle listing.
 */
template <bool is_weighted, typename NodeAddedCallbackType>
class LocalDegreeDirectedGraph {
protected:
    // local to global (input graph) id mapping
    std::vector<node> local_to_global_id;
    // map from input graph to local id
    std::unordered_map<node, node> global_to_local_id;

    // for the local graph: outgoing neighbors
    std::vector<node> head;
    // for the local graph: weight of outgoing edges, only used for weighted graphs (otherwise
    // empty)
    std::vector<double> head_weight;

    struct local_node {
        /** first index in head where neighbors of the node are stored */
        std::size_t first_head;
        /** index after the last index in head where neighbors of the node are stored */
        std::size_t last_head;
    };

    // stores for every local node the information where the outgoing neighbors begins and ends
    std::vector<local_node> head_info;
    // stores the degree of every local node in the local graph
    std::vector<count> degree;

    std::vector<node> current_local_neighbor_ids;

    const NetworKit::Graph &G;

    // data structures only for neighbors of a node
    std::vector<bool> is_neighbor;
    std::vector<double> neighbor_weight;

    // callback to call whenever a node is added
    NodeAddedCallbackType node_added_callback;

public:
    /**
     * Initialize the local graph as empty graph.
     *
     * @param G The graph from which nodes can be added
     * @param node_added_callback Callback that is called whenever a node is added
     */
    LocalDegreeDirectedGraph(const NetworKit::Graph &G, NodeAddedCallbackType node_added_callback)
        : G(G), node_added_callback(node_added_callback) {}

    ~LocalDegreeDirectedGraph() = default;

    /**
     * Map a local node id to a global one
     */
    node toGlobal(node lu) const { return local_to_global_id[lu]; }

    /**
     * Ensure that the given global node id exists in the local graph
     *
     * @param u The global node id to check
     * @return The local node id
     */
    node ensureNodeExists(node u) {
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

            // Iterate over all neighbors in the local graph
            for (index i = my_begin; i < my_end;) {
                node ln = head[i];
                edgeweight weight = is_weighted ? head_weight[i] : 1;

                auto &n_degree = degree[ln];
                auto &n_info = head_info[ln];

                // Before into the neighbors of ln, check if it can get rid of any
                // neighbors due to being a higher-degree node now.
                for (index ni = n_info.first_head; ni < n_info.last_head;) {
                    node nn = head[ni];
                    if (n_degree <= degree[nn]) {
                        ++ni;
                    } else {
                        // Move edge to node nn
                        // Store last neighbor of ln at position ni instead, i.e., remove nn
                        head[ni] = head[--n_info.last_head];

                        if (is_weighted) {
                            // Copy weight to neighbor list of nn
                            head_weight[head_info[nn].last_head] = head_weight[ni];
                            // Copy weight of last neighbor of nn at position ni
                            head_weight[ni] = head_weight[n_info.last_head];
                        }

                        // Store ln as neighbor at the last position of the neighbor list of nn
                        head[head_info[nn].last_head++] = ln;
                    }
                }

                if (my_degree < n_degree) {
                    ++i;
                } else {
                    // Add local_id to ln as neighbor instead of adding ln as neighbor here
                    if (is_weighted) {
                        head_weight[n_info.last_head] = weight;
                    }

                    head[n_info.last_head++] = local_id;

                    // This removes one neighbor
                    --my_end;
                    // Move the last neighbor at the current position
                    head[i] = head[my_end];
                    if (is_weighted) {
                        head_weight[i] = head_weight[my_end];
                    }

                    assert(ln + 1 == local_id
                           || head_info[ln].last_head <= head_info[ln + 1].first_head);
                }
            }

            head_info.emplace_back(local_node{my_begin, my_end});

            is_neighbor.push_back(false);

            if (is_weighted)
                neighbor_weight.push_back(0);

            node_added_callback(u, local_id, weightedDegree);

            return local_id;
        } else {
            return gtl_it->second;
        }
    }

    /**
     * Iterate over all triangles that include the global node @a u
     *
     * @param u The node from which triangles shall be listed
     * @param callback The callback to call. Must take two node ids of the triangle as well as three
     * edge weights as parameters.
     */
    template <typename F>
    void forTrianglesOf(node u, F callback) {
        forTrianglesOf(
                       u, [](node, edgeweight) {}, callback);
    }

    /**
     * Iterate over all triangles that include the global node @a u
     *
     * @param u The node from which triangles shall be listed
     * @param initNeighborsCallback The callback to call for every local neighbor (local id) as they are collected
     * @param triangleCallback The callback to call for each
     * triangle. Must take two node ids of the triangle as well as
     * three edge weights as parameters.
     */
    template <typename Fi, typename Ft>
    void forTrianglesOf(node u, Fi initNeighborsCallback, Ft triangleCallback) {
        if (G.degree(u) == 0)
            return;

        current_local_neighbor_ids.clear();
        current_local_neighbor_ids.reserve(G.degree(u));

        G.forNeighborsOf(u, [&](node, node v, edgeweight weight) {
            node local_id = ensureNodeExists(v);
            current_local_neighbor_ids.emplace_back(local_id);
            is_neighbor[local_id] = true;

            if (is_weighted) {
                neighbor_weight[local_id] = weight;
            }

            initNeighborsCallback(local_id, weight);
        });

        for (node lv : current_local_neighbor_ids) {
            for (index ti = head_info[lv].first_head; ti < head_info[lv].last_head; ++ti) {
                node y = head[ti];

                if (is_neighbor[y]) {
                    edgeweight weight_uv = 1.0;
                    edgeweight weight_uy = 1.0;
                    edgeweight weight_vy = 1.0;
                    if (is_weighted) {
                        weight_uv = neighbor_weight[lv];
                        weight_uy = neighbor_weight[y];
                        weight_vy = head_weight[ti];
                    }

                    triangleCallback(lv, y, weight_uv, weight_uy, weight_vy);
                }
            }
        }

        for (node v : current_local_neighbor_ids) {
            is_neighbor[v] = false;
        }
    }

    /**
     * Iterate over local neighbors of the node previously passed to forTrianglesOf
     *
     * @param callback The call back that shall be called for each neighbor, needs to take a node
     * and an edgeweight.
     */
    template <typename F>
    void forLocalNeighbors(F callback) {
        for (node v : current_local_neighbor_ids) {
            if (is_weighted) {
                callback(v, neighbor_weight[v]);
            } else {
                callback(v, 1.0);
            }
        }
    }
};

} // namespace NetworKit
#endif
