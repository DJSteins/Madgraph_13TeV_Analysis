
#include "HEPUtils/Event.h"
#include "TMVA/Reader.h"

#ifdef DANIELS_COMP
#include "RestFrames/RestFrames.hh"
#endif

using namespace std;

namespace Gambit {
  namespace ColliderBit {
    class Analysis_ATLAS_13TeV_1LEPStop_36invfb {
    private:
      // defines an enumerator that has corresponding strings
#define SRMAP(X) \
        X(tN_med) \
        X(tN_high) \
        X(bWN) \
        X(bC2x_diag) \
        X(bC2x_med) \
        X(bCbv) \
        X(bffN) \
        X(DM_low_loose) \
        X(DM_low) \
        X(DM_high) \
        X(tN_diag_low)  \
        X(tN_diag_med)  \
        X(tN_diag_high) \
        X(bCsoft_diag)  \
        X(bCsoft_med) \
        X(bCsoft_high)
#define f(x) x,
#define g(x) #x,
      enum SignalRegion {SRMAP(f)};
      const std::vector<std::string> SignalRegionStrings = {SRMAP(g)};
#undef g
#undef f
#undef SRMAP

      // stores the number of times that an event was detected for a signal region
      map<string, int> signalRegionCounts;

      vector<double> BDTLowScores;
      vector<double> BDTMedScores;
      vector<double> BDTHighScores;

      // a map to store the values of variables in so we can save them in a CSV and plot them later
      map<string, vector<double>> varResults;
      vector<string> varNames = {
        "DMTAlpha",
        "tN_diag_low_MET",
        "tN_diag_low_mT",
        "topHadronicMass",
        "topLeptonic200Mass",
        "dPhiLepNu200",
        "dPhiRjetLep",
        "preRjetPt",
        "preDPhiRjetLep",
        "preBLepMass",
        "preTopHadronic",
        "preMTopLep200",
        "preDMTAlpha",
        "preDPhiLepNu200",
        "preMet",
        "preMt",
        "preNjets",
        "preNBjets",

        "mTopChi2",
        "tN_diag_med_MET",
        "tN_diag_med_mT",
        "tN_diag_med_htSigMiss",
        "tN_diag_med_dRbl",
        "tN_diag_med_dPhiTTbar",
        "dPhiHadTopMet",
        "tN_diag_med_nJets",
        "jet_pt[2]",
        "jet_pt[3]",
        "tN_diag_med_preMTopChi2",

        "tN_diag_high_mT",
        "tN_diag_high_dPhiTTbar",
        "tN_diag_high_dRbl",
        "jetPt3",
        "jetPt4",
        "tN_diag_high_mTopChi2",
        "njv",
        "mTS",
        "rISR",
        "dPhiISR",
        "tN_diag_high_nJets",
        "tN_diag_high_MET"
      };

      int cutFlowIndex;
      vector<int> cutFlowVector;
      vector<string> cutFlowVector_str;
      unsigned NCUTS;

      vector<HEPUtils::Particle*> baselineLeptons;
      vector<HEPUtils::Particle*> signalLeptons;
      vector<HEPUtils::Jet*> signalJets;
      vector<double> signalJER;
      HEPUtils::P4 metVec;
      size_t nJets;
      size_t nBJets;
      double Met;
      double absDPhiJMet0;
      double absDPhiJMet1;
      double mT2Tau;
      double jetPt1;
      double jetPt2;
      double jetPt3;
      double jetPt4;
      double jetPt5;
      double MetPerp;
      double HtSigMiss;
      double mT;
      double amT2;
      double dRbl;
      double topReclM;
      double dPhiMetLep;
      double minDPhiMetBJet;
      double pTLepOverMet;
      double bJetPt1;
      double bJetPt2;
      double WReclM;
      double Wpt;
      double dRbb;
      double mW;
      double mTop;
      bool preselSoftLep;
      bool jet1Btag;

      // RJR variables
      int njv;
      double rISR;
      double dphiISRI;
      double mTS;

      string BDTWeightsLocation;
      string BDTWeightsFileNum;

      unique_ptr<RestFrames::LabRecoFrame> LAB;
      unique_ptr<RestFrames::DecayRecoFrame> CM;
      unique_ptr<RestFrames::DecayRecoFrame> S;
      unique_ptr<RestFrames::VisibleRecoFrame> ISR;
      unique_ptr<RestFrames::VisibleRecoFrame> V;
      unique_ptr<RestFrames::InvisibleRecoFrame> I;
      unique_ptr<RestFrames::InvisibleGroup> INV;
      unique_ptr<RestFrames::CombinatoricGroup> VIS;
      unique_ptr<RestFrames::SetMassInvJigsaw> InvMass;
      unique_ptr<RestFrames::MinMassesCombJigsaw> SplitVis;

      /**
        * Updates the cutFlowVector and cutFlowVector_str using the given vectors
        *
        * @param checkVector the requirements for the signal region
        * @param baseLabel the name of the signal region
        * @param checkStrVector the descriptions of each of the cuts
        */
      void updateCutFlowVector(
        const vector<bool>& checkVector,
        const string& baseLabel,
        const vector<string>& checkStrVector);


      /**
       * Fills the given signal vectors with the appropriate particles.
       *
       * @param baselineParticles the particles to search through
       * @param softSignalParticles filled with all particles in baselineParticles
       * @param signalParticles filled with all particles in baselineParticles
       * @param softSignalLeptons filled with particles from baselineParticles with pT > pTLimit
       * @param pTLimit the max allowed pT for the signal particles/leptons
       */
      void fillSignalParticles(
        vector<HEPUtils::Particle*> *baselineParticles,
        vector<HEPUtils::Particle*> *softSignalParticles,
        vector<HEPUtils::Particle*> *signalParticles,
        vector<HEPUtils::Particle*> *softSignalLeptons,
        double pTLimit);

      /**
       * Selects the two indices for the bJets based on signal jet btags
       *
       * @return a vector of {bjet1 index, bjet2 index}
       */
      vector<int> bJetSelection();

      /**
       * Saves the cutflow table to a file.
       */
      void saveCutFlow();

      /**
       * Adds SignalRegionData for a given signal region's information. The info can be found in tables 22-25 in the
       * paper.
       *
       * @param sr_label the label for the signal region
       * @param n_observed the number of observed events
       * @param n_background the number of background events
       * @param background_sys the background systematic uncertainty
       * @param signal_sys the signal systematic uncertainty
       */
      void addSignalRegionData(
        const SignalRegion& sr_label,
        int n_observed,
        double n_background,
        double background_sys,
        double signal_sys);

      /**
       * Saves BDT scores to a file with the given fileName
       *
       * @param vec the scores to save
       * @param fileName the name of the output file
       */
      void save1DVector(vector<double> vec, const string &fileName, int numEvents);

      /**
       * Saves the simulation cross-section to a txt file
       */
      void saveXSec();

      /**
       * Sets a lot of the analysis variables to zero that are set as class variables
       */
      void initVars();

      /**
       * Performs a chi squared selection to pick the most hadronic top quark and best b-lepton
       *
       * @param topHadronic the hadronic top quark to set
       * @param bLepton the b-lepton to set
       * @param mostBjetLike a vector of the most b-like jets
       * @param signalNotBjetLike a vector of the most not b-like signal jets
       */
      void chi2Selection(HEPUtils::P4* topHadronic, HEPUtils::P4* bLepton, vector<HEPUtils::Jet*> mostBjetLike, vector<HEPUtils::Jet*> signalNotBjetLike);

      /**
       * Gets the jet comb
       *
       * @param bJetNums
       * @param mostBjetLike
       * @return a vector of jet comb indices
       */
      vector<int> getJetComb(vector<int> bJetNums, const vector<HEPUtils::Jet*>& mostBjetLike, const HEPUtils::BinnedFn2D<double>& _resJets2D);

      /**
       * Fills the signalNotBjet vectors
       * @param signalNotBjet a vector to fill
       * @param signalNotBjetLike a vector to fill
       * @param bJetNums
       */
      void getSignalNotBjets(vector<HEPUtils::Jet*> *signalNotBjet, vector<HEPUtils::Jet*> *signalNotBjetLike, vector<int> bJetNums);

      /**
       * Calculates HtSigMiss
       * @return the value of HtSigMiss
       */
      double calcHtSigMiss();

      /**
       * Searches for the best absDPhiJiMet from the four possibilities given
       * @param absDPhiJMet an array of possibilities
       * @return absDPhiJiMet
       */
      double getAbsDPhiJiMet(const double absDPhiJMet[4]);

      /**
       * Does recursive jigsaw operations to get the RJR variables
       */
      void doRJRCalculations();

      /**
       * Assigns some of the commonly used momenta and tags from the signal jets
       * @param signalBJets a vector of signal b-jets
       */
      void assignJetValues(vector<HEPUtils::Jet*> signalBJets);

      /**
       * Calculates amT2
       * @param baseLepton the first lepton in baselineLeptons
       * @param mostBjetLike a vector of most b-jet like jets
       * @return the value of amT2
       */
      double calcAmt2(HEPUtils::Particle* baseLepton, vector<HEPUtils::Jet*> mostBjetLike);

      /**
       * Calculates MetPerp
       * @param topChi2 a chi-squared reconstructed top
       * @param top1 another top
       * @return the value of MetPerp
       */
      double calcMetPerp(const HEPUtils::P4& topChi2, const HEPUtils::P4& top1);

      /**
       * Calculates mT2Tau
       * @param baselineTaus a vector of baseline tau particles
       * @param baseLepton the first lepton in baselineLeptons
       * @return the value of mT2Tau
       */
      double calcMT2Tau(vector<HEPUtils::Particle*> baselineTaus, HEPUtils::Particle* baseLepton);

    public:
      /**
       * The constructor that should initialize some variables and the Recursive Jigsaw Lab frames
       */
      Analysis_ATLAS_13TeV_1LEPStop_36invfb();

      /**
       * Performs the main part of the analysis
       * @param event an event contain particle and jet information
       */
      void analyze(const HEPUtils::Event* event) ;

      /**
       * Stores results into SignalRegionData and adds it to _results
       */
      void collect_results(int numEvents) ;

      /**
       * Assigns the values for absDPhiJMet based on the number of signal jets we have
       * @param absDPhiJMet an array to fill with the values
       */
      void assignAbsDPhiJMet(double absDPhiJMet[4]);

      /**
       * Increments the given signalRegionName's count by 1 in the signalRegionCounts map based on the given boolean
       * @param signalRegion the key corresponding to a signal region in signalRegionCounts
       * @param isAccepted the boolean that determines whether the signal region is accepted
       */
      void acceptSignalRegion(SignalRegion signalRegion, bool isAccepted);

      /**
       * Adds a variable to the given reader by converting the given double to a new float(value)
       * @param reader the reader to update
       * @param varName the name of the variable in the xml weights file
       * @param value the value of the variable
       */
      inline void addBDTVariable(TMVA::Reader* reader, const string& varName, double value);

      /**
       * Loads the BDT Weights file number from BDTWeightsLocation + "BDT_File_Number.txt"
       */
      void loadBDTFileNumber();
      /**
       * Resets some variables
       */
      void clear() ;

    protected:

      /**
       * Performs an overlap removal procedure for two leptons as specified in the paper in table 3
       * @param lep1vec the first vector of leptons
       * @param lep2vec the second vector of leptons
       * @param DeltaRMax the maximum delta R
       */
      void LeptonLeptonOverlapRemoval(vector<HEPUtils::Particle*> *lep1vec, vector<HEPUtils::Particle*> *lep2vec, double DeltaRMax);

      /**
       * Performs an overlap removal procedure where object 1 is a lepton and object 2 is a jet as specified in table 3
       * @param jetvec a vector of jets
       * @param lepvec a vector of leptons
       * @param DeltaRMax the maximum delta R
       */
      void JetLeptonOverlapRemoval(vector<HEPUtils::Jet*> *jetvec, vector<HEPUtils::Particle*> *lepvec, double DeltaRMax);

      /**
       * Performs and overlap removal procedure where object 1 is a jet and object 2 is a lepton as specified in table 3
       * @param lepvec a vector of leptons
       * @param jetvec a vector of jets
       */
      void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*> *lepvec, vector<HEPUtils::Jet*> *jetvec);

      /**
       * Loops through each of the pT values from the signal jets and checks them against the requirements from the
       * paper that this file is analyzing.
       *
       * @param pTRequirements the list of pT requirements specified in the paper
       * @param signalJetsPts the signal jet pt values
       */
      bool checkJetPt(vector<int> pTRequirements, vector<double> signalJetPts);

      /**
       * Gets the baseline leptons that from the given particles that satisfy the given requirements.
       *
       * @param eventParticles the specific particles to sort through
       * @param minPt the minimum pT requirement
       * @param maxAbsEta the maximum absolute eta requirement
       * @return a vector of baseline particles
       */
      vector<HEPUtils::Particle*> filterPtEta(vector<HEPUtils::Particle*> eventParticles, double minPt, double maxAbsEta);

      /**
       * Gets the jets that satisfy the given requirements and stores them in a new vector.
       *
       * @param typedJets the specific signal jets we want to look at
       * @param minPt the minimum pT requirement
       * @param maxAbsEta the maximum absolute eta requirement
       * @param bTagSetting the btag to store if the requirements are satisfied
       * @return a vector of signal jets
       */
      vector<HEPUtils::Jet*> getSignalTypeJets(vector<HEPUtils::Jet*> typedJets, double minPt, double maxAbsEta, bool bTagSetting);

      /**
       * Chooses b-jets and non b-jets based on a tagging efficiency, pT, and absolute eta
       * @param event the event to analyze
       * @param bJets a vector to store b-jets in
       * @param nonBJets a vector to store non b-jets in
       */
      void getBJets(const HEPUtils::Event* event, vector<HEPUtils::Jet*>* bJets, vector<HEPUtils::Jet*>* nonBJets);

    };
  }
}