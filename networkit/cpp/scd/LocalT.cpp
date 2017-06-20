#include "LocalT.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <limits>

namespace NetworKit {

LocalT::LocalT(const Graph& G) : SelectiveCommunityDetector(G) {
}

namespace {

template <typename NodeAddedCallbackType>
class LocalGraph {
private:
	// local to global (input graph) id mapping
	std::vector<node> local_to_global_id;
	// map from input graph to local id
	std::unordered_map<node, node> global_to_local_id;

	// for the local graph: outgoing neighbors
	std::vector<node> head;

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

	// callback to call whenever a node is added
	NodeAddedCallbackType node_added_callback;
public:
	LocalGraph(const NetworKit::Graph &G, NodeAddedCallbackType node_added_callback) : G(G), node_added_callback(node_added_callback) {
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

			size_t my_degree = 0;

			index my_begin = head.size();
			index my_end = my_begin;
			index nh = my_end; // next head if all potential out neighbors were inserted

			G.forEdgesOf(u, [&](node v) {
				auto it = global_to_local_id.find(v);

				if (it != global_to_local_id.end()) {
					++my_degree;

					auto& n_degree = degree[it->second];
					++n_degree;

					head.push_back(it->second);
					++my_end;
				}

				++nh;
			});

			assert(my_end == head.size());
			assert(nh >= my_end);

			head.resize(nh);

			degree.push_back(my_degree);

			for (index i = my_begin; i < my_end;) {
				node ln = head[i];

				auto& n_degree = degree[ln];
				auto& n_info = head_info[ln];

				// before pushing into n_edges check if we can get rid of any neighbors.
				for (index ni = n_info.first_head; ni < n_info.last_head;) {
					node nn = head[ni];
					if (n_degree <= degree[nn]) {
						++ni;
					} else {
						head[ni] = head[--n_info.last_head];
						head[head_info[nn].last_head++] = ln;
					}
				}

				if (my_degree < n_degree) {
					++i;
				} else {
					head[n_info.last_head++] = local_id;

					--my_end;
					head[i] = head[my_end];
					//assert(head_info[ln].last_head <= head_info[ln+1].first_head);
				}
			}

			head_info.emplace_back(local_node {my_begin, my_end});

			//head_info.back().last_head = head.size();

			is_neighbor.push_back(false);
			node_added_callback(u, local_id);

			return local_id;
		} else {
			return gtl_it->second;
		}
	}

	template <typename F>
	void forTrianglesOf(node u, F callback) {
		if (G.degree(u) == 0)
			return;

		current_local_neighbor_ids.clear();
		current_local_neighbor_ids.reserve(G.degree(u));

		G.forNeighborsOf(u, [&](node v) {
			node local_id = ensureNodeExists(v);
			current_local_neighbor_ids.emplace_back(local_id);
			is_neighbor[local_id] = true;
		});


		for (node lv : current_local_neighbor_ids) {
			for (index ti = head_info[lv].first_head; ti < head_info[lv].last_head; ++ti) {
				node y = head[ti];

				if (is_neighbor[y]) {
					callback(lv, y);
				}

			}
		}

		for (node v : current_local_neighbor_ids) {
			is_neighbor[v] = false;
		}
	}
};


}

std::set<node> LocalT::expandOneCommunity(const std::set<node>& s) {
	// global data structures
	// result community
	std::set<node> result;

	// This algorithm creates a local graph. It contains the community and its shell.
	// stores the sum of the similarities of all nodes to nodes in the community
	std::vector<double> node_internal_triangles;
	std::vector<double> node_external_triangles;
	std::vector<double> node_semi_internal_triangles;;

	// indicates for every local node if it is in the community (true) or in the shell (false)
	std::vector<bool> in_result;

	std::vector<bool> in_shell;

	count internal_triangles = 0, external_triangles = 0;

	auto addNode = [&](node, node) {
		node_internal_triangles.push_back(0);
		node_external_triangles.push_back(0);
		node_semi_internal_triangles.push_back(0);
		in_result.push_back(false);
		in_shell.push_back(false);
	};

	LocalGraph<decltype(addNode)> local_graph(G, addNode);

	std::unordered_set<node> shell;


	auto updateShell = [&] (node u) {
		local_graph.forTrianglesOf(u, [&](node lv, node lw) {
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

		G.forNeighborsOf(u, [&](node v) {
			node s = local_graph.ensureNodeExists(v);

			if (!in_shell[s] && !in_result[s]) {
				shell.insert(s);
				in_shell[s] = true;

				local_graph.forTrianglesOf(local_graph.toGlobal(s), [&](node lv, node lw) {
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
			local_graph.forTrianglesOf(u, [&](node lv, node lw) {
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
		assert(debug_internal_triangles/3 == internal_triangles);

		for (node ls : shell) {
			count s_external = 0, s_internal = 0, s_semi_internal = 0;
			local_graph.forTrianglesOf(local_graph.toGlobal(ls), [&](node lv, node lw) {
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
	for (node u: s) {
		node lu = local_graph.ensureNodeExists(u);
		result.insert(u);
		in_result[lu] = true;
		updateShell(u);
	}

	auto get_score = [](count int_triangles, count ext_triangles) -> count {
		return std::max<int64_t>(0, int_triangles * (static_cast<int64_t>(int_triangles) - static_cast<int64_t>(ext_triangles)));
	};

	// expand (main loop)
	node uMax = none;
	do {

		uMax = none;
		count best_score = get_score(internal_triangles, external_triangles);
		count best_external_triangles = none;

		for (node lv : shell) {
			count new_internal_triangles = internal_triangles + node_internal_triangles[lv];
			count new_external_triangles = external_triangles + node_external_triangles[lv] - node_semi_internal_triangles[lv];

			count new_score = get_score(new_internal_triangles, new_external_triangles);

			if (new_score > best_score || (new_score == best_score && new_external_triangles < best_external_triangles)) {
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
