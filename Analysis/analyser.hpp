
/* analysis module as abstract base class */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

template<class T>
class Analyser {
  public:
    
    Analyser() {};
    virtual ~Analyser() {};
    virtual void perform_analysis() = 0;
    virtual void write_analysis_results() = 0;
  
};
