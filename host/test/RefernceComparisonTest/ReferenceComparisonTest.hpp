#pragma once

#include <memory>
#include <vector>
#include "../Ssv.hpp"
#include "../../types/HavacHit.hpp"

using std::shared_ptr;
using std::vector;
using std::make_shared;

void compareHardwareSoftwareHits(vector<HavacHit> hardwareHits, shared_ptr<vector<ReferenceSsvHit>> softwareHits);