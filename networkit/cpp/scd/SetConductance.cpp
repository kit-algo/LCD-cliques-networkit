#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/scd/SetConductance.hpp>

namespace NetworKit {
SetConductance::SetConductance(const Graph &G, const std::set<node> &community)
    : G(&G), community(&community) {}

void SetConductance::run() {
    hasRun = false;

    Aux::SignalHandler handler;

    double cut = 0;
    double com_volume = 0;

    for (node u : *community) {
        if (G->hasNode(u)) {
            G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                if (community->count(v) == 0) {
                    cut += ew;
                }
                com_volume += ew;
            });
        }
    }

    double total_volume = G->totalEdgeWeight() * 2;

    double rest_volume = total_volume - com_volume;

    conductance = 1.0;

    if (com_volume > 0 && rest_volume > 0) {
        conductance = cut / std::min(com_volume, rest_volume);
    }

    hasRun = true;
}

double SetConductance::getConductance() const {
    assureFinished();
    return conductance;
}

} // namespace NetworKit
