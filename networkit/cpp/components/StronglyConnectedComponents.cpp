/*
 * StronglyConnectedComponents.cpp
 *
 *  Created on: 01.06.2014
 *      Authors: Klara Reichard <klara.reichard@gmail.com>
 *               Marvin Ritter <marvin.ritter@gmail.com>
 *               Obada Mahdi <omahdi@gmail.com>
 *               Eugenio Angriman <angrimae@hu-berlin.de>
 */

// networkit-format

#include <stack>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

StronglyConnectedComponents::StronglyConnectedComponents(const Graph &G) : G(&G) {
    if (!G.isDirected())
        WARN("The input graph is undirected, use ConnectedComponents for more efficiency.");
}

void StronglyConnectedComponents::run() {
    const auto n = G->upperNodeIdBound();

    // The component of every node is initially undefined.
    component.clear();
    component.resize(n, none);

    // Stack used to determine the strongly connected components from the root nodes.
    std::vector<node> stack;
    stack.reserve(n);
    std::vector<unsigned char> onStack(n, 0);

    // Depth of a node in the dfs
    std::vector<count> depth(n, none);
    // At index `u`: smallest id among reachable nodes.
    std::vector<node> lowLink(n, none);
    // Last visited child during the dfs.
    std::vector<node> lastVisited(n, none);

    count curDepth = 0, visitedNodes = 0;
    componentIndex = 0;

    std::stack<std::pair<node, Graph::NeighborIterator>> dfsStack;

    // Set the depth of node v and push it onto the stacks.
    auto visit = [&](const node v) -> void {
        depth[v] = curDepth;
        lowLink[v] = curDepth;
        ++curDepth;
        stack.emplace_back(v);
        onStack[v] = 1;
        dfsStack.push({v, G->neighborRange(v).begin()});
    };

    auto strongConnect = [&](const node u) -> void {
        visit(u);

        do {
            const auto v = dfsStack.top().first;
            if (lastVisited[v] != none) {
                lowLink[v] = std::min(lowLink[v], lowLink[lastVisited[v]]);
                lastVisited[v] = none;
            }

            // Iter points to the first neighbor of v, or to the last visited dfs child of v.
            auto &iter = dfsStack.top().second;

            // Iterate over the neighbors of v from either the first, or the last visited one.
            for (; iter != G->neighborRange(v).end(); ++iter) {
                const auto w = *iter;

                // w not visited yet, visit it and continue the exploration from w.
                if (depth[w] == none) {
                    visit(w);
                    lastVisited[v] = w;
                    break;
                }
                // w already visited, if it is on the stack is part of the same component.
                if (onStack[w] && depth[w] < lowLink[v]) {
                    lowLink[v] = depth[w];
                }
            }

            // Check if all neighbors of v have been visited.
            if (iter == G->neighborRange(v).end()) {
                // All neighbors of v have been visited, pop v.
                dfsStack.pop();

                // v is a root node, generate new component.
                if (lowLink[v] == depth[v]) {

                    const auto stackSize = stack.size();
                    node w = none;
                    do {
                        w = stack.back();
                        stack.pop_back();
                        onStack[w] = 0;
                        component[w] = componentIndex;
                    } while (w != v);

                    ++componentIndex;
                    visitedNodes += (stackSize - stack.size());
                }
            }
        } while (!dfsStack.empty());
    };

    for (node u = 0; u < n; ++u) {
        if (depth[u] != none)
            continue;
        strongConnect(u);
        // Check if all nodes have been assigned to a component.
        if (visitedNodes == G->numberOfNodes())
            break;
    }

    hasRun = true;
}

} // namespace NetworKit
