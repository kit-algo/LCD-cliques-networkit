#include "LocalTightnessExpansion.h"

#include <algorithm>
#include <unordered_map>
#include <limits>
#include "../auxiliary/DynamicMinMaxHeap.h"

using Aux::MaxHeap;

namespace NetworKit {

LocalTightnessExpansion::LocalTightnessExpansion(const Graph& G, double alpha) : SelectiveCommunityDetector(G), alpha(alpha) {
}

namespace {

#ifndef NDEBUG
/*
 * @deprecated Only used for debugging.
 */
double weightedEdgeScore(const Graph&G, node u, node v, edgeweight ew, double uDegree, const std::unordered_map<node, double>& uNeighbors) {
	double vDegree = 1.0; // w(u, u) = 1 per their definition
	double nom = ew; // w(v, v) * w(v, u)

	G.forNeighborsOf(v, [&] (node, node w, edgeweight vwWeight) {
		vDegree += vwWeight * vwWeight;
		auto uw = uNeighbors.find(w);
		if (uw != uNeighbors.end())
			nom += vwWeight * uw->second;
		else if (w == u)
			nom += vwWeight; // w(u, u) * w(u, v)
	});

	vDegree = std::sqrt(vDegree);

	double denom = vDegree * uDegree;

	return nom /denom;
}
#endif

template <bool is_weighted, typename NodeAddedCallbackType>
class LocalGraph {
private:
	// local to global (input graph) id mapping
	std::vector<node> local_to_global_id;
	// map from input graph to local id
	std::unordered_map<node, node> global_to_local_id;

	// for the local graph: outgoing neighbors
	std::vector<node> head;
	// for the local graph: weight of outgoing edges, only used for weighted graphs (otherwise empty)
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
	std::vector<double> triangle_sum;
	std::vector<double> neighbor_weight;

	// callback to call whenever a node is added
	NodeAddedCallbackType node_added_callback;
public:
	LocalGraph(const NetworKit::Graph &G, NodeAddedCallbackType node_added_callback) : G(G), node_added_callback(node_added_callback) {
		if (G.isWeighted() != is_weighted) {
			throw std::runtime_error("Error, weighted/unweighted status of input graph does not match is_weighted template parameter");
		}
	}

	node toGlobal(node lu) const {
		return local_to_global_id[lu];
	}

	node ensureNodeExists(node u) {
		auto gtl_it = global_to_local_id.find(u);

		if (gtl_it == global_to_local_id.end()) {
			node local_id = local_to_global_id.size();
			local_to_global_id.push_back(u);
			global_to_local_id[u] = local_id;

			double weightedDegree = 1.0; // w(u, u) = 1
			//count degU = G.degree(u);

			size_t my_degree = 0;

			index my_begin = head.size();
			index my_end = my_begin;
			index nh = my_end; // next head if all potential out neighbors were inserted

			G.forEdgesOf(u, [&](node, node v, edgeweight weight) {
				auto it = global_to_local_id.find(v);

				if (it != global_to_local_id.end()) {
					++my_degree;

					auto& n_degree = degree[it->second];
					++n_degree;

					head.push_back(it->second);
					if (is_weighted) {
						head_weight.push_back(weight);
					}

					++my_end;
				}

				weightedDegree += weight * weight;
				++nh;
			});

			weightedDegree = std::sqrt(weightedDegree);

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

				auto& n_degree = degree[ln];
				auto& n_info = head_info[ln];

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

					//assert(head_info[ln].last_head <= head_info[ln+1].first_head);
				}
			}

			head_info.emplace_back(local_node {my_begin, my_end});

			//head_info.back().last_head = head.size();

			is_neighbor.push_back(false);
			triangle_sum.push_back(0);

			if (is_weighted) {
				neighbor_weight.push_back(0);
			}

			node_added_callback(u, local_id);

			return local_id;
		} else {
			return gtl_it->second;
		}
	}

	template <typename F>
	void forNeighborsWithTrianglesOf(node u, F callback) {
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
			triangle_sum[local_id] = 2 * weight;
		});


		for (node lv : current_local_neighbor_ids) {
			// count (weighted) triangles using only out-neighbors
			double triangles_lv = 0;
			for (index ti = head_info[lv].first_head; ti < head_info[lv].last_head; ++ti) {
				node y = head[ti];

				if (is_neighbor[y]) {
					if (is_weighted) {
						triangle_sum[y] += neighbor_weight[lv] * head_weight[ti];
						triangles_lv += neighbor_weight[y] * head_weight[ti];
					} else {
						triangle_sum[y] += 1;
						triangles_lv += 1;
					}
				}

			}

			triangle_sum[lv] += triangles_lv;
		}

		for (node v : current_local_neighbor_ids) {
			double w = 1;
			if (is_weighted) w = neighbor_weight[v];
			callback(v, w, triangle_sum[v]);
			is_neighbor[v] = false;
		}
	}
};



template <bool is_weighted>
std::set<node> expandSeedSet_internal(const Graph&G, const std::set<node>& s, double alpha) {
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
	MaxHeap<double> shell;

 // 	auto globalIsOutNeighbor = [&](node u, count degU, node v) {
 // 		count degV = G.degree(v);
 // 		return (degV > degU);// || (degV == degU && v > u);
 // 	};

	double internal_similarity = 0;
	double external_similarity = 0;

	auto addNode = [&](node u, node) {
		double wd = 1;
		if (is_weighted) {
			G.forNeighborsOf(u, [&](node, node, edgeweight w) {
				wd += w*w;
			});
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


	auto updateShell = [&] (node u, node lu, bool noverify = false) {
#ifndef NDEBUG
		std::unordered_map<node, double> uNeighbors;
		G.forNeighborsOf(u, [&](node, node v, edgeweight ew) {
				uNeighbors.insert(std::make_pair(v, ew));
		});

		if (in_shell[lu]) {
			double debug_u_external_similarity = .0;
			double debug_u_internal_similarity = .0;
			local_graph.forNeighborsWithTrianglesOf(u, [&](node lv, edgeweight, double triangle_sum) {
				double score_uv = 0.0;
				double denom = weighted_degree[lv] * weighted_degree[lu];
				if (denom > 0.0)  {
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

		local_graph.forNeighborsWithTrianglesOf(u, [&](node lv, edgeweight weight, double triangle_sum) {
		// collect counts and compute scores
			double denom = weighted_degree[lv] * weighted_degree[lu];
			double score_uv = triangle_sum / denom;

#ifndef NDEBUG
			{
				double debug_score = weightedEdgeScore(G, u, local_graph.toGlobal(lv), weight, weighted_degree[lu], uNeighbors);
				assert(score_uv == debug_score);
			}

			if (shell.contains(lv)) {
				assert(shell.getKey(lv) == node_internal_similarity[lv]);
			}
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

				shell.pushOrIncreaseKey(lv, node_internal_similarity[lv]);

				if (!in_shell[lv]) {
					in_shell[lv] = true;
					new_shell_nodes.push_back(lv);
				} else {
					node_external_similarity[lv] -= score_uv;
				}

#ifndef NDEBUG
				if (false && !noverify) {
					// FIXME adapt this sanity check for new measure
					double debugWeight = 0;
					G.forNeighborsOf(local_graph.toGlobal(lv), [&](node, node y, edgeweight ew) {
						if (result.count(y)) {
							debugWeight += ew;
						}
					});
					assert(debugWeight == node_internal_similarity[lv]);
				}
#endif
			}
		});

		for (node s : new_shell_nodes) {
			local_graph.forNeighborsWithTrianglesOf(local_graph.toGlobal(s), [&](node lv, edgeweight, double triangle_sum) {
				if (!in_result[lv]) {
					double score_uv = 0.0;
					double denom = weighted_degree[lv] * weighted_degree[s];
					if (denom > 0.0)  {
						score_uv = triangle_sum / denom;
					}

					node_external_similarity[s] += score_uv;
				}
			});
		}
	};

	// init community with seed set
	for (node u: s) {
		node lu = local_graph.ensureNodeExists(u);
		result.insert(u);
		in_result[lu] = true;
		updateShell(u, lu);
	}

	// expand (main loop)
	while (!shell.empty()) {
		node uMax = shell.pop();
		node gUMax = local_graph.toGlobal(uMax);
#ifndef NDEBUG
		double internal_similarity_new = internal_similarity + 2 * node_internal_similarity[uMax];
		double external_similarity_new = external_similarity - node_internal_similarity[uMax] + node_external_similarity[uMax];

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
		if (external_similarity / internal_similarity - (alpha * node_external_similarity[uMax] - node_internal_similarity[uMax]) / (2 * node_internal_similarity[uMax]) > 0) {
			result.insert(gUMax);
			in_result[uMax] = true;
			updateShell(gUMax, uMax);

			assert(std::abs(internal_similarity - internal_similarity_new) < 0.000001);
			assert(std::abs(external_similarity - external_similarity_new) < 0.000001);
		}
	}
	return result;
}

} // internal namespace

std::set<node> LocalTightnessExpansion::expandOneCommunity(const std::set<node>& s) {
	if (G.isWeighted()) {
		return expandSeedSet_internal<true>(G, s, alpha);
	} else {
		return expandSeedSet_internal<false>(G, s, alpha);
	}
};


} /* namespace NetworKit */
