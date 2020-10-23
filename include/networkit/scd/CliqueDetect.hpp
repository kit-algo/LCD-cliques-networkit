#ifndef NETWORKIT_SCD_CLIQUE_DETECT_HPP_
#define NETWORKIT_SCD_CLIQUE_DETECT_HPP_

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

/**
 * The CliqueDetect algorithm. First finds the largest clique in the seed node's
 * neighborhood.
 *
 * The algorithm can handle weighted graphs.
 */
class CliqueDetect : public NetworKit::SelectiveCommunityDetector {

public:
    /**
     * Construct a Cliquedetect object.
     *
     * @param[in] G The graph to detect communities on
     */
    CliqueDetect(const Graph &G);

    /**
     * Expands a single seed node/vertex into a maximal clique.
     *
     * @param[in] s the seed node
     * @return A community of the seed node
     */
    std::set<node> expandOneCommunity(node seed) override;

    /**
     * Detect a single clique for the given seed nodes.
     *
     * The resulting community is a clique iff the seeds form a clique.
     * Otherwise, only the added nodes form a clique that is fully connected
     * to the seed nodes.
     *
     * @param seeds The seeds for the community.
     * @return The found community as set of nodes.
     */
    std::set<node> expandOneCommunity(const std::set<node> &seeds) override;

protected:
    template <typename F>
    std::vector<node> getMaximumWeightClique(const std::vector<node> &nodes, F is_seed) const;

    Graph createSubgraphFromNodes(const std::vector<node> &nodes) const;
};

} // namespace NetworKit

#endif // NETWORKIT_SCD_CLIQUE_DETECT_HPP_
