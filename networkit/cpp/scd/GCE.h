/* GCE.h
 *
 * Created on: 06.05.2013
 * Author: cls
 */


#ifndef GCE_H_
#define GCE_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"
#include "../auxiliary/SetIntersector.h"


namespace NetworKit {


/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
class GCE: public NetworKit::SelectiveCommunityDetector {

public:

	GCE(const Graph& G, std::string objective);


	/**
	 * @param[in]	seeds	seed nodes
	 *
	 * @param[out]		community as a set of nodes
	 */
	std::set<node> expandOneCommunity(const std::set<node>& s) override;

	// inherit method from parent class.
	using SelectiveCommunityDetector::expandOneCommunity;

protected:

    std::string objective;    // name of objective function


};

} /* namespace NetworKit */
#endif
