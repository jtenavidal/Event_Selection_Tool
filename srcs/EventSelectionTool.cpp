#include "../include/EventSelectionTool.h"
#include <iostream>

namespace selection{
  
  void LoadEventList(const std::string file_name){
  
    TFile f(file_name.c_str());

    TTree *event_tree = (TTree*) f.Get("event_tree");

  
  }

} // namespace: selection

