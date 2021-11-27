#include <iostream>
#include "highway_cover_labelling.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc > 1) {
    int k = atoi(argv[3]);

    cout << "Loading Graph..." << std::endl;
    HighwayLabelling hl(argv[2], k);

    int topk[k];
    hl.SelectLandmarks_HD(topk);

    if(string(argv[1]).compare("construct_labelling") == 0) {
      cout << "Constructing Highway Cover Labelling..." << std::endl;
      hl.BuildIndex(topk);
      
      hl.StoreLabelling(argv[4]); //storing labelling to disk
      hl.deallocate(0);
    } else if(string(argv[1]).compare("update_labelling") == 0) {
      hl.LoadLabelling_Full(argv[4], topk, k); //loading labelling from disk

      cout << "Updating Highway Cover Labelling..." << std::endl;
      hl.UpdateLabelling(argv[5], atoi(argv[6]));
      
      hl.StoreLabelling(argv[4]); //storing labelling to disk after update
      hl.deallocate(0);
    } else if (string(argv[1]).compare("query_labelling") == 0) {
      hl.LoadLabelling_Pruned(argv[4]); //loading labelling from disk
      hl.RemoveLandmarks(topk);
      
      cout << "Performing Distance Queries..." << std::endl;
      hl.QueryDistance(argv[5], argv[6]);
      
      hl.deallocate(1);
    }
  }

  return 0;
}
