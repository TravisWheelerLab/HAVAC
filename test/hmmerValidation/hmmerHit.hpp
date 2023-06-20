#pragma once

#include <string>
#include <cstdint>
#include <vector>
#include <sstream>

using std::string;
using std::vector;

struct HmmerHitLine {
  string target;
  string accession;
  string query;
  uint32_t hmmfrom;
  uint32_t hmmto;
  uint32_t alifrom;
  uint32_t alito;
  uint32_t envfrom;
  uint32_t envto;
  uint32_t seqlen;
  string strand;
  float eValue;
  float score;
  float bias;

  string toString() {
    std::stringstream ss;
    ss << "target: " << target << ", acc: " << accession << ", query: " << query << ", hmmfrom: " << hmmfrom <<
      ", hmmto: " << hmmto << ", envfrom: " << envfrom << ", envto: " << envto << ", eval: " << eValue;

    return ss.str();
  }
};



HmmerHitLine getHitFromLine(string lineBuffer);
vector<HmmerHitLine> getHitsFromFile(string fileSrc);