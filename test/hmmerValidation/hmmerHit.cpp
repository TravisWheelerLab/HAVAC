#include "hmmerHit.hpp"
#include <stdexcept>
#include <sstream>
#include <fstream>


HmmerHitLine getHitFromLine(string lineBuffer) {
  if (lineBuffer[0] == '#') {
    throw std::logic_error("line buffer given was a comment line (starts with #)");
  }

  HmmerHitLine hit;
  std::istringstream stream(lineBuffer);
  string dummyString;
  stream >> hit.target;
  stream >> dummyString; //skip the sequence accession!
  stream >> hit.query >> hit.accession >> hit.hmmfrom >> hit.hmmto >> hit.alifrom >>
    hit.alifrom >> hit.envfrom >> hit.envto >> hit.seqlen >> hit.strand >> hit.eValue >> hit.score >> hit.bias;

  return hit;
}


vector<HmmerHitLine> getHitsFromFile(string fileSrc) {
  vector<HmmerHitLine> hmmerHits;
  string lineBuffer;
  std::fstream fstream;
  fstream.open(fileSrc, std::ios::in);
  //read the first two (dummy) lines
  std::getline(fstream, lineBuffer);
  std::getline(fstream, lineBuffer);

  //read the hits, stopping on the comments at the end.
  while (true) {
    std::getline(fstream, lineBuffer);
    if (lineBuffer.size() == 0 || lineBuffer[0] == '#') {
      break;
    }

    HmmerHitLine hit = getHitFromLine(lineBuffer);
    hmmerHits.push_back(hit);
  }

  return hmmerHits;
}