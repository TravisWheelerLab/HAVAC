#pragma once

#include <memory>
#include <vector>
#include "Ssv.hpp"
#include "../../types/VerifiedHit.hpp"

using std::shared_ptr;
using std::vector;
using std::make_shared;

void compareHardwareSoftwareHits(shared_ptr<vector<VerifiedHit>> hardwareHits, shared_ptr<vector<ReferenceSsvHit>> softwareHits);