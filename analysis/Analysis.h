#include "NPVUserAnalysis.h"
// #include "MyDetector.h"
namespace user_analysis {
  class Analysis : public nptool::VUserAnalysis {
   public:
    Analysis(){};
    ~Analysis(){};

   public:
    void Init();
    void TreatEvent();
    // if this method return false, the event is discarded
    // Caution : this change the size of the tree,
    // so the friend mecanism no longer work
    bool FillOutputCondition() { return true; };
    void End();

   private:
    //  std::shared_ptr<my_detector::MyDetector> mydetector;
  };
} // namespace user_analysis
