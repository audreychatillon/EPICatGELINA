#include "NPVUserAnalysis.h"
#include "EpicDetector.h"

namespace user_analysis {
  class Analysis : public nptool::VUserAnalysis {
   public:
    Analysis(){};
    ~Analysis(){};

   public:
    void Init();
    void TreatEvent();
    void SetDataOutput(std::shared_ptr<nptool::VDataOutput>);
    // if this method return false, the event is discarded
    // Caution : this change the size of the tree,
    // so the friend mecanism no longer work
    bool FillOutputCondition() { return true; };
    void End();

   private:
    std::shared_ptr<epic::EpicDetector> epic;
  };
} // namespace user_analysis
