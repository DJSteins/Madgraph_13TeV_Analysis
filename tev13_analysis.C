#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <random>
#include <iostream>
#include <fstream>
#include "gambit/ColliderBit/analyses/AnalysisUtil.hpp"
#include "tev13_analysis.h"

#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "TRandom3.h"

//HEPUtils
#include "HEPUtils/FastJet.h"
#include "HEPUtils/Jet.h"
#include "HEPUtils/Event.h"
#include "HEPUtils/Particle.h"
#include "HEPUtils/BinnedFn.h"
#include "HEPUtils/MathUtils.h"
#include "MCUtils/PIDCodes.h"

//Other files
#include "gambit/ColliderBit/lester_mt2_bisect.h"
#include "gambit/ColliderBit/mt2_bisect.h"
#include <cmath>

//Delphes and root stuff
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif


using namespace std;

/* The ATLAS 1 lepton direct stop analysis

   Based on: https://arxiv.org/abs/1711.11520

   Code by Martin White (based on ATLAS public code snippet on HepData)

   Modified and extended by Daniel Steiner

   KNOWN ISSUES

   1) Have not added the BDT signal regions (despite having BDT code from ATLAS). They cover a specific kinematic region where the m_stop - m_chi1 mass difference is m_top, which we already know Pythia does badly with.

   2) We have no equivalent of the ATLAS fakeJER method. Am assuming a 3% JER on every jet for now.

   3) Have used TLorentzVectors for boosting. Could probably be done without ROOT?

*/

double Gambit::Random::draw()
{
  Utils::specialised_threadsafe_rng<mt19937_64> ultralocal_rng;
  Utils::threadsafe_rng* local_rng = &ultralocal_rng;
  return (*local_rng)();
}

namespace Gambit {
  namespace ColliderBit {
// Need two different functions here for use with sort
    bool sortByPT_1l(HEPUtils::Jet* jet1, HEPUtils::Jet* jet2)
    {
      return (jet1->pT() > jet2->pT());
    }

    bool sortByPT_1l_Particle(HEPUtils::Particle* p1, HEPUtils::Particle* p2)
    {
      return (p1->pT() > p2->pT());
    }

    bool sortByPT_1l_sharedptr(const shared_ptr<HEPUtils::Jet>& jet1, const shared_ptr<HEPUtils::Jet>& jet2)
    {
      return sortByPT_1l(jet1.get(), jet2.get());
    }

    double calcMT_1l(const HEPUtils::P4& jetMom, const HEPUtils::P4& metMom)
    {
      double met = sqrt(metMom.px() * metMom.px() + metMom.py() * metMom.py());
      double dphi = metMom.deltaPhi(jetMom);
      double mt = sqrt(2 * jetMom.pT() * met * (1 - cos(dphi)));
      return mt;
    }

    double sqr(double value)
    {
      return value * value;
    }

/**
 * Gets deltaphi of two phis and returns it. no absolute value here
 * @param a
 * @param b
 * @return
 */
    double dPhi(double a, double b)
    {
      double rtn = a - b;
      rtn = fmod(rtn, 2*M_PI);
      assert(rtn >= -2*M_PI && rtn <= 2*M_PI);
      if (rtn == 0) return 0;
      if (rtn > M_PI) rtn -= 2*M_PI;
      if (rtn <= -M_PI) rtn += 2*M_PI;
      assert(rtn > -M_PI && rtn <= M_PI);
      return rtn;
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::loadBDTFileNumber()
    {
      ifstream stream;
      stream.open(BDTWeightsLocation + "BDT_File_Number.txt");
      string key;
      string value;
      getline(stream, key, ':');
      getline(stream, value, '\n');
      // should be a 1 character number
      if (value.length() == 1)
      {
        BDTWeightsFileNum = value;
        cout << "Using BDT File Number: " << BDTWeightsFileNum << endl;
        stream.close();
      }
      else
      {
        cout << "Invalid value found in ATLAS_13TeV_1LEPStop_36invfb /BDT_File_Number.txt";
        stream.close();
        exit(0);
      }
    }

    Analysis_ATLAS_13TeV_1LEPStop_36invfb::Analysis_ATLAS_13TeV_1LEPStop_36invfb()
    {
      for (const string& key : SignalRegionStrings)
      {
        signalRegionCounts.insert(make_pair(key, 0));
      }

      for (const string& varName : varNames)
      {
        varResults.insert(make_pair(varName, vector<double>({})));
      }

      NCUTS = 150;

      cutFlowVector = vector<int>(NCUTS, 0);
      cutFlowVector_str = vector<string>(NCUTS, "");

      BDTWeightsLocation = "/home/dsteiner/gambit/ColliderBit/data/ATLAS_13TeV_1LEPStop_36invfb/";
      loadBDTFileNumber();

      // RestFrames initialisation

      LAB = make_unique<RestFrames::LabRecoFrame>("LAB", "lab");
      CM = make_unique<RestFrames::DecayRecoFrame>("CM", "cm");
      S = make_unique<RestFrames::DecayRecoFrame>("S", "s");
      ISR = make_unique<RestFrames::VisibleRecoFrame>("ISR", "isr");
      V = make_unique<RestFrames::VisibleRecoFrame>("V", "v");
      I = make_unique<RestFrames::InvisibleRecoFrame>("I", "i");

      // Connect the frames
      LAB->SetChildFrame(*CM);
      CM->AddChildFrame(*ISR);
      CM->AddChildFrame(*S);
      S->AddChildFrame(*V);
      S->AddChildFrame(*I);

      // Initialize the tree
      LAB->InitializeTree();

      // Define groups
      INV = make_unique<RestFrames::InvisibleGroup>("INV", "inv");
      INV->AddFrame(*I);
      VIS = make_unique<RestFrames::CombinatoricGroup>("VIS", "vis");
      VIS->AddFrame(*ISR);
      VIS->SetNElementsForFrame(*ISR, 1, false);
      VIS->AddFrame(*V);
      VIS->SetNElementsForFrame(*V, 0, false);

      // set the invisible system mass to zero
      InvMass = make_unique<RestFrames::SetMassInvJigsaw>("InvMass", "kSetMass");
      INV->AddJigsaw(*InvMass);

      // define the rule for partitioning objects between "ISR" and "V"
      SplitVis = make_unique<RestFrames::MinMassesCombJigsaw>("CombPPJigsaw", "kMinMasses");
      VIS->AddJigsaw(*SplitVis);
      // "0" group (ISR)
      SplitVis->AddFrame(*ISR, 0);
      // "1" group (V + I)
      SplitVis->AddFrame(*V, 1);
      SplitVis->AddFrame(*I, 1);

      LAB->InitializeAnalysis();
    }

    struct ClusteringHistory : public FJNS::PseudoJet::UserInfoBase
    {
      enum Status
      {
        GOOD,
        JET_TOO_SMALL,
        JET_TOO_LARGE,
        TOO_MANY_ITERATIONS,
        NONE,
      };

      struct Step
      {
        double pt;
        double r;
        size_t constit;
        Status status;
      };

      size_t id;  // a per-event unique jet id that is needed for the event dump
      vector<Step> steps;

      static ClusteringHistory* AddStep(ClusteringHistory& history, const Step& step)
      {
        auto newHistory = new ClusteringHistory(history);
        newHistory->steps.push_back(step);
        return newHistory;
      }
    };

// Return the history of a PseudoJet object, handling all the ugly casting.
    ClusteringHistory& GetHistory(const FJNS::PseudoJet& jet)
    {
      return *dynamic_cast<ClusteringHistory*>(jet.user_info_shared_ptr().get());
    }

    static vector<FJNS::PseudoJet> SortedByNConstit(vector<FJNS::PseudoJet> jets)
    {
      sort(jets.begin(), jets.end(),
           [](const FJNS::PseudoJet& a, const FJNS::PseudoJet& b) {
             if (a.constituents().size() != b.constituents().size()) {
               return a.constituents().size() > b.constituents().size();
             }
             return a.pt() > b.pt();
           });

      return jets;
    }

    inline double optimalRadius(const double pT, const double m) { return 2 * m / pT; }

    inline double minRadius(const double pT, const double m) { return optimalRadius(pT, m) - 0.3; }

    inline double maxRadius(const double pT, const double m) { return optimalRadius(pT, m) + 0.5; }

    inline double chi2(const double observation, const double expectation)
    {
      return pow(observation - expectation, 2) / expectation;
    }

    pair<bool, FJNS::PseudoJet> RecursiveRecluster(const FJNS::PseudoJet& candidate, double candRadius,
                                                   const double mass, size_t step)
    {
      if (minRadius(candidate.pt(), mass) > candRadius)
      {
        GetHistory(candidate).steps.back().status = ClusteringHistory::JET_TOO_SMALL;
        return make_pair(false, candidate);
      }
      else if (maxRadius(candidate.pt(), mass) < candRadius)
      {
        const double newR = max(maxRadius(candidate.pt(), mass), candRadius / 2.);
        GetHistory(candidate).steps.back().status = ClusteringHistory::JET_TOO_LARGE;

        if (step > 10)
        {
          GetHistory(candidate).steps.back().status = ClusteringHistory::TOO_MANY_ITERATIONS;
          return make_pair(false, candidate);
        }

        FJNS::JetDefinition jetDef(FJNS::antikt_algorithm, newR);
        auto cs = new FJNS::ClusterSequence(candidate.constituents(), jetDef);

        vector<FJNS::PseudoJet> reclusteredJets;
        reclusteredJets = SortedByNConstit(cs->inclusive_jets());

        if (reclusteredJets.empty())
        {
          delete cs;
          return make_pair(false, FJNS::PseudoJet());
        }

        cs->delete_self_when_unused();
        auto newCandidate = reclusteredJets[0];

        auto newHistory = ClusteringHistory::AddStep(
          GetHistory(candidate),
          {newCandidate.pt(), newR, newCandidate.constituents().size(), ClusteringHistory::NONE});
        newCandidate.set_user_info(newHistory);

        return RecursiveRecluster(newCandidate, newR, mass, step + 1);
      }
      else
      {
        GetHistory(candidate).steps.back().status = ClusteringHistory::GOOD;
        return make_pair(true, candidate);
      }
    }


    HEPUtils::P4 reclusteredParticle(vector<HEPUtils::Jet*> jets, vector<HEPUtils::Jet*> bjets,
                           const double mass, const bool useBJets)
    {
      HEPUtils::P4 p;
      double r0 = 3.0;

      vector<HEPUtils::Jet*> usejets;
      for(HEPUtils::Jet* jet : jets)
      {
        usejets.push_back(jet);
      }

      if (useBJets && !bjets.empty())
      {
        for(HEPUtils::Jet* bjet : bjets)
        {
          usejets.push_back(bjet);
        }
      }

      vector<FJNS::PseudoJet> initialJets;

      for (HEPUtils::Jet* jet : usejets)
      {
        FJNS::PseudoJet Pjet(jet->mom().px(), jet->mom().py(), jet->mom().pz(), jet->mom().E());
        initialJets.push_back(Pjet);
      }

      FJNS::JetDefinition jetDef(FJNS::antikt_algorithm, r0);
      FJNS::ClusterSequence cs(initialJets, jetDef);

      auto candidates = FJNS::sorted_by_pt(cs.inclusive_jets());

      vector<FJNS::PseudoJet> selectedJets;
      selectedJets.reserve(candidates.size());
      vector<FJNS::PseudoJet> badJets;
      badJets.reserve(candidates.size());

      size_t i = 0;
      for (auto& cand : candidates)
      {
        auto history = new ClusteringHistory();
        history->id = i;
        history->steps.push_back({cand.pt(), r0, cand.constituents().size(), ClusteringHistory::NONE});
        cand.set_user_info(history);
        ++i;
      }

      for (const auto& cand : candidates)
      {
        bool selected = false;
        FJNS::PseudoJet jet;

        tie(selected, jet) = RecursiveRecluster(cand, r0, mass, 0);

        if (selected)
          selectedJets.push_back(jet);
        else
          badJets.push_back(jet);
      }

      if (selectedJets.empty())
      {
        return p;
      }

      vector<shared_ptr<HEPUtils::Jet>> aoSelectedJets;
      for (const FJNS::PseudoJet& j : selectedJets) aoSelectedJets.push_back(make_shared<HEPUtils::Jet>(HEPUtils::mk_p4(j)));

      sort(aoSelectedJets.begin(), aoSelectedJets.end(), sortByPT_1l_sharedptr);
      p = aoSelectedJets[0]->mom();

      return p;
    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::LeptonLeptonOverlapRemoval(
      vector<HEPUtils::Particle*> *lep1vec,
      vector<HEPUtils::Particle*> *lep2vec,
      double DeltaRMax)
    {

      //Routine to do jet-lepton check
      //Discards jets if they are within DeltaRMax of a lepton

      vector<HEPUtils::Particle*> Survivors;

      for (HEPUtils::Particle* lep1 : *lep1vec)
      {
        bool overlap = false;
        for (HEPUtils::Particle* lep2 : *lep2vec)
        {
          double dR = lep1->mom().deltaR_eta(lep2->mom());

          if(fabs(dR) <= DeltaRMax) overlap = true;
        }
        if (overlap) continue;
        Survivors.push_back(lep1);
      }
      *lep1vec = Survivors;
    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::JetLeptonOverlapRemoval(
      vector<HEPUtils::Jet*> *jetvec,
      vector<HEPUtils::Particle*> *lepvec,
      double DeltaRMax)
    {
      //Routine to do jet-lepton check
      //Discards jets if they are within DeltaRMax of a lepton

      vector<HEPUtils::Jet*> Survivors;

      for (HEPUtils::Jet* jet : *jetvec)
      {
        bool overlap = false;
        for (HEPUtils::Particle* lep : *lepvec)
        {
          double dR = jet->mom().deltaR_eta(lep->mom());

          if(fabs(dR) <= DeltaRMax) overlap = true;
        }
        if(overlap) continue;
        Survivors.push_back(jet);
      }
      *jetvec = Survivors;
    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::LeptonJetOverlapRemoval(vector<HEPUtils::Particle*> *lepvec, vector<HEPUtils::Jet*> *jetvec)
    {
      //Routine to do lepton-jet check
      //Discards leptons if they are within dR of a jet as defined in analysis paper

      vector<HEPUtils::Particle*> Survivors;

      for (HEPUtils::Particle* lep : *lepvec)
      {
        bool overlap = false;
        for (HEPUtils::Jet* jet : *jetvec)
        {
          double dR = jet->mom().deltaR_eta(lep->mom());
          double DeltaRMax = max(0.1, min(0.4, 0.04 + 10 / lep->mom().pT()));

          if(fabs(dR) <= DeltaRMax) overlap = true;
        }
        if(overlap) continue;
        Survivors.push_back(lep);
      }
      *lepvec = Survivors;
    }


    bool Analysis_ATLAS_13TeV_1LEPStop_36invfb::checkJetPt(vector<int> pTRequirements, vector<double> signalJetPts)
    {
      if (pTRequirements.size() > signalJetPts.size())
      {
        return false;
      }
      int i = 0;
      for (int pT : pTRequirements)
      {
        if (signalJetPts[i++] <= pT)
        {
          // end here if any of the requirements are not met
          return false;
        }
      }
      // all requirements met
      return true;
    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::updateCutFlowVector(
      const vector<bool>& checkVector,
      const string& baseLabel,
      const vector<string>& checkStrVector)
    {
      unsigned i;
      // when a requirement isn't satisfied, this needs to be set to false
      int cutFlowIncrement = 1;
      // loop through the checks and increment cuts while the requirements are satisfied
      for (i = 0; i < checkVector.size(); i++)
      {
        cutFlowVector_str[cutFlowIndex] = baseLabel + checkStrVector.at(i);
        if (checkVector.at(i)) {
          cutFlowVector[cutFlowIndex] += cutFlowIncrement;
        } else {
          cutFlowIncrement = 0;
        }
        cutFlowIndex++;
      }
    }


    vector<HEPUtils::Particle*> Analysis_ATLAS_13TeV_1LEPStop_36invfb::filterPtEta(
      vector<HEPUtils::Particle *> eventParticles,
      double minPt,
      double maxAbsEta)
    {
      vector<HEPUtils::Particle*> outParticles;
      for (HEPUtils::Particle* particle : eventParticles)
      {
        if (particle->pT() > minPt && particle->abseta() < maxAbsEta)
        {
          outParticles.push_back(particle);
        }
      }
      return outParticles;
    }


    vector<HEPUtils::Jet*> Analysis_ATLAS_13TeV_1LEPStop_36invfb::getSignalTypeJets(
      vector<HEPUtils::Jet*> typedJets,
      double minPt,
      double maxAbsEta,
      bool bTagSetting)
    {
      vector<HEPUtils::Jet*> signalTypeJets;
      for (HEPUtils::Jet* jetType : typedJets)
      {
        if(jetType->pT() > minPt && fabs(jetType->eta()) < maxAbsEta)
        {
          jetType->set_btag(bTagSetting);
          signalJets.push_back(jetType);
          signalTypeJets.push_back(jetType);
        }
      }
      return signalTypeJets;
    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::fillSignalParticles(
      vector<HEPUtils::Particle*> *baselineParticles,
      vector<HEPUtils::Particle*> *softSignalParticles,
      vector<HEPUtils::Particle*> *signalParticles,
      vector<HEPUtils::Particle*> *softSignalLeptons,
      double pTLimit)
    {
      for (HEPUtils::Particle* particle : *baselineParticles)
      {
        softSignalParticles->push_back(particle);
        softSignalLeptons->push_back(particle);
        if(particle->pT() > pTLimit)
        {
          signalParticles->push_back(particle);
          signalLeptons.push_back(particle);
        }
      }
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::getSignalNotBjets(
      vector<HEPUtils::Jet*> *signalNotBjet,
      vector<HEPUtils::Jet*> *signalNotBjetLike,
      vector<int> bJetNums)
    {
      int i;
      for (i = 0; i < (int) signalJets.size(); ++i)
      {
        if (!signalJets[i]->btag())
        {
          signalNotBjet->push_back(signalJets[i]);
        }
        if (i == bJetNums[0] || i == bJetNums[1])
          continue;
        signalNotBjetLike->push_back(signalJets[i]);
      }
    }

    double Analysis_ATLAS_13TeV_1LEPStop_36invfb::calcHtSigMiss()
    {
      /* calculate vecHtMiss */
      HEPUtils::P4 vecHtMiss, leptonHtMiss;
      for (HEPUtils::Particle *baselineLepton : baselineLeptons) {
        vecHtMiss -= baselineLepton->mom();
        leptonHtMiss -= baselineLepton->mom();
      }

      double Ht = 0.0;
      for (HEPUtils::Jet *jet : signalJets) {
        vecHtMiss -= jet->mom();
        Ht += jet->pT();
      }

      TRandom3 myRandom;
      myRandom.SetSeed((ULong_t) signalJets[0]->pT());

      int PEs = 100;
      double ETmissmean = 0, ETmissRMS = 0;
      size_t i;
      for (int j = 0; j < PEs; ++j) {
        double jetHtx = leptonHtMiss.px();
        double jetHty = leptonHtMiss.py();

        for (i = 0; i < signalJets.size(); ++i) {
          jetHtx -= myRandom.Gaus(signalJets[i]->mom().px(), signalJets[i]->mom().px() * signalJER[i]);
          jetHty -= myRandom.Gaus(signalJets[i]->mom().py(), signalJets[i]->mom().px() * signalJER[i]);
        }
        double ETtemp = sqrt(jetHtx * jetHtx + jetHty * jetHty);
        ETmissmean += ETtemp;
        ETmissRMS += ETtemp * ETtemp;
      }

      ETmissmean = ETmissmean / PEs;
      double sigmaAbsHtMiss = sqrt((ETmissRMS / PEs) - ETmissmean * ETmissmean);
      return (ETmissmean - 100.) / sigmaAbsHtMiss;
    }

    double Analysis_ATLAS_13TeV_1LEPStop_36invfb::getAbsDPhiJiMet(const double absDPhiJMet[4])
    {
      double absDPhiJiMet = absDPhiJMet[0];
      size_t i;
      for (i = 1; i < 4; i++)
      {
        if (absDPhiJMet[i] < absDPhiJiMet)
        {
          absDPhiJiMet = absDPhiJMet[i];
        }
      }
      return absDPhiJiMet;
    }

    vector<int> Analysis_ATLAS_13TeV_1LEPStop_36invfb::bJetSelection()
    {
      // create containers with exactly 2 jets being considered to be b-jets and the inverse
      int bJet1 = -1, bJet2 = -1;
      unsigned i;
      for (i = 0; i < nJets; ++i)
      {
        if (!signalJets[i]->btag()) continue;
        if (bJet1 == -1)
          bJet1 = i;
        else if (bJet2 == -1)
        {
          bJet2 = i;
          break;
        }
      }
      if (bJet2 == -1)
      {
        for (i = 0; i < nJets; ++i)
        {
          if (signalJets[i]->btag()) continue;
          if (bJet1 == -1)
            bJet1 = i;
          else if (bJet2 == -1)
          {
            bJet2 = i;
            break;
          }
        }
      }
      return {bJet1, bJet2};
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::getBJets(const HEPUtils::Event* event, vector<HEPUtils::Jet*>* bJets, vector<HEPUtils::Jet*>* nonBJets)
    {
      /// @note We assume that b jets have previously been 100% tagged
      const vector<double> a = {0, 10.};
      const vector<double> b = {0, 10000.};
      const vector<double> c = {0.77}; // set b-tag efficiency to 77%
      HEPUtils::BinnedFn2D<double> _eff2d(a,b,c);
      for (HEPUtils::Jet* jet : event->jets())
      {
        bool hasTag = has_tag(_eff2d, jet->eta(), jet->pT());
        if (jet->pT() > 20. && jet->abseta() < 4.9)
        {
          if(jet->btag() && hasTag && jet->abseta() < 2.5)
          {
            bJets->push_back(jet);
          }
          else
          {
            nonBJets->push_back(jet);
          }
        }
      }
    }

    vector<int> Analysis_ATLAS_13TeV_1LEPStop_36invfb::getJetComb(
      vector<int> bJetNums,
      const vector<HEPUtils::Jet*>& mostBjetLike,
      const HEPUtils::BinnedFn2D<double>& _resJets2D)
    {
      vector<double> signalBJER;
      for (HEPUtils::Jet *jet : mostBjetLike) {
        signalBJER.push_back(_resJets2D.get_at(jet->abseta(), jet->pT()));
      }

      auto chi2min = DBL_MAX;
      double f;
      vector<int> jetComb = {0, 0, 0};
      int i, j, k;
      for (i = 0; i < (int) signalJets.size(); i++)
      {
        if (i == bJetNums[0] || i == bJetNums[1])
          continue;
        for (j = i + 1; j < (int) signalJets.size(); j++)
        {
          if (j == bJetNums[0] || j == bJetNums[1])
            continue;
          for (k = 0; k < (int) mostBjetLike.size() && k < 2; k++) {
            HEPUtils::P4 signalJetSum = signalJets[i]->mom() + signalJets[j]->mom();
            double particleMass1 = (signalJetSum + mostBjetLike[k]->mom()).m();
            double particleMass2 = signalJetSum.m();
            double sumOfSquares2 = pow(signalJER[i], 2) + pow(signalJER[j], 2);
            double sumOfSquares1 = sumOfSquares2 + pow(signalBJER[k], 2);
            f = pow(particleMass1 - mTop, 2) / (pow(particleMass1, 2) * sumOfSquares1) +
                pow(particleMass2 - mW, 2) / (pow(particleMass2, 2) * sumOfSquares2);
            if (f < chi2min) {
              chi2min = f;
              jetComb[0] = i;
              jetComb[1] = j;
              jetComb[2] = k;
            }
          }
        }
      }
      return jetComb;
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::chi2Selection(
      HEPUtils::P4* topHadronic,
      HEPUtils::P4* bLepton,
      vector<HEPUtils::Jet*> mostBjetLike,
      vector<HEPUtils::Jet*> signalNotBjetLike)
    {
      // run a chi squared minimiztaion to pick the best hadronic top and leptonic b quarks
//  unsigned i, j;
      double chi2_1, chi2_2;
      auto chi2min = DBL_MAX;
      const double mW = 80.385;
      const double mTop = 173.21;
      for (unsigned int i = 1; i < signalNotBjetLike.size(); ++i) {
        //1-jet W
        chi2_1 =
          pow(signalNotBjetLike[i]->mass() - mW, 2) / mW +
          pow((*signalNotBjetLike[i] + *mostBjetLike[0]).m() - mTop, 2) / mTop;
        chi2_2 =
          pow(signalNotBjetLike[i]->mass() - mW, 2) / mW +
          pow((*signalNotBjetLike[i] + *mostBjetLike[1]).m() - mTop, 2) / mTop;

        if (chi2_1 < chi2_2 && chi2_1 < chi2min) {
          chi2min = chi2_1;
          *topHadronic = *signalNotBjetLike[i] + *mostBjetLike[0];
          *bLepton = *mostBjetLike[1];
        } else if (chi2_2 < chi2_1 && chi2_2 < chi2min) {
          chi2min = chi2_2;
          *topHadronic = *signalNotBjetLike[i] + *mostBjetLike[1];
          *bLepton = *mostBjetLike[0];
        }
        for (unsigned int j = i+1; j < signalNotBjetLike.size(); ++j) {

          //2-jet W
          chi2_1 = pow((*signalNotBjetLike[i] + *signalNotBjetLike[j]).m() - mW, 2) / mW +
                   pow((*signalNotBjetLike[i] + *signalNotBjetLike[j] + *mostBjetLike[0]).m() - mTop, 2) / mTop;
          chi2_2 = pow((*signalNotBjetLike[i] + *signalNotBjetLike[j]).m() - mW, 2) / mW +
                   pow((*signalNotBjetLike[i] + *signalNotBjetLike[j] + *mostBjetLike[1]).m() - mTop, 2) / mTop;

          if (chi2_1 < chi2_2 && chi2_1 < chi2min) {
            chi2min = chi2_1;
            *topHadronic = *signalNotBjetLike[i] + *signalNotBjetLike[j] + *mostBjetLike[0];
            *bLepton = *mostBjetLike[1];
          } else if (chi2_2 < chi2_1 && chi2_2 < chi2min) {
            chi2min = chi2_2;
            *topHadronic = *signalNotBjetLike[i] + *signalNotBjetLike[j] + *mostBjetLike[1];
            *bLepton = *mostBjetLike[0];
          }
        }
      }
//  for (i = 1; i < signalNotBjetLike.size(); i++)
//  {
//    //1-jet W
//    HEPUtils::P4 firstMostBjetLike = *mostBjetLike[0];
//    HEPUtils::P4 secondMostBjetLike = *mostBjetLike[1];
//    HEPUtils::P4 thisSignalNotBjetLike = *signalNotBjetLike[i];
//    HEPUtils::P4 firstJetSum = thisSignalNotBjetLike + firstMostBjetLike;
//    HEPUtils::P4 secondJetSum = thisSignalNotBjetLike + secondMostBjetLike;
//    double jetChi2 = chi2(thisSignalNotBjetLike.m(), mW);
//    chi2_1 = jetChi2 + chi2(firstJetSum.m(), mTop);
//    chi2_2 = jetChi2 + chi2(secondJetSum.m(), mTop);
//
//    if (chi2_1 < chi2_2 && chi2_1 < chi2min) {
//      chi2min = chi2_1;
//      *topHadronic = firstJetSum;
//      *bLepton = secondMostBjetLike;
//    } else if (chi2_2 < chi2_1 && chi2_2 < chi2min) {
//      chi2min = chi2_2;
//      *topHadronic = secondJetSum;
//      *bLepton = firstMostBjetLike;
//    }
//    for (j = i + 1; j < signalNotBjetLike.size(); j++)
//    {
//      //2-jet W
//      HEPUtils::P4 subSignalNotBjetLike = *signalNotBjetLike[j];
//      HEPUtils::P4 subJetSum = thisSignalNotBjetLike + subSignalNotBjetLike;
//      HEPUtils::P4 firstSubJetSum = subJetSum + firstMostBjetLike;
//      HEPUtils::P4 secondSubJetSum = subJetSum + secondMostBjetLike;
//      double subJetChi2 = chi2(subJetSum.m(), mW);
//      chi2_1 = subJetChi2 + chi2(firstSubJetSum.m(), mTop);
//      chi2_2 = subJetChi2 + chi2(secondSubJetSum.m(), mTop);
//
//      if (chi2_1 < chi2_2 && chi2_1 < chi2min) {
//        chi2min = chi2_1;
//        *topHadronic = firstSubJetSum;
//        *bLepton = secondMostBjetLike;
//      } else if (chi2_2 < chi2_1 && chi2_2 < chi2min) {
//        chi2min = chi2_2;
//        *topHadronic = secondSubJetSum;
//        *bLepton = firstMostBjetLike;
//      }
//    }
//  }
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::doRJRCalculations()
    {
      // Do RestFrames Calculations
      LAB->ClearEvent();
      vector<RestFrames::RFKey> jetID;

      for (HEPUtils::Jet* jet : signalJets)
      {
        TLorentzVector fourVector;
        fourVector.SetPtEtaPhiM(jet->pT(), 0.0, jet->phi(), jet->mass());
        jetID.push_back(VIS->AddLabFrameFourVector(fourVector));
      }

      TVector3 ETMiss;
      ETMiss.SetXYZ(metVec.px(), metVec.py(), 0.0);
      INV->SetLabFrameThreeVector(ETMiss);

      bool analysisSuccessful = LAB->AnalyzeEvent();
      if (!analysisSuccessful)
      {
        cout << "Something went wrong..." << endl;
      }

      for (const RestFrames::RFKey& jetIDKey : jetID)
      {
        if (VIS->GetFrame(jetIDKey).IsSame(*V))
        { // sparticle group
          njv++;
        }
      }
      // need at least one jet associated with MET-side of event
      if (njv > 0)
      {
        TVector3 vP_ISR = ISR->GetFourVector(*CM).Vect();
        TVector3 vP_I   = I->GetFourVector(*CM).Vect();

        rISR = fabs(vP_I.Dot(vP_ISR.Unit())) / vP_ISR.Mag();
        dphiISRI = fabs(vP_ISR.DeltaPhi(vP_I));
        mTS = S->GetMass();
      }
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::assignJetValues(vector<HEPUtils::Jet*> signalBJets)
    {
      jetPt1 = signalJets[0]->pT();
      jetPt2 = signalJets[1]->pT();
      nJets > 2 ? jetPt3 = signalJets[2]->pT() : 0.0;
      nJets > 3 ? jetPt4 = signalJets[3]->pT() : 0.0;
      nJets > 4 ? jetPt5 = signalJets[4]->pT() : 0.0;
      jet1Btag = signalJets[0]->btag();

      nBJets > 0 ? bJetPt1 = signalBJets[0]->pT() : 0.0;
      nBJets > 1 ? bJetPt2 = signalBJets[1]->pT() : 0.0;
    }

    double Analysis_ATLAS_13TeV_1LEPStop_36invfb::calcAmt2(HEPUtils::Particle* baseLepton, vector<HEPUtils::Jet*> mostBjetLike)
    {
      HEPUtils::P4 lepton_plus_bjet0 = baseLepton->mom() + mostBjetLike[0]->mom();
      HEPUtils::P4 lepton_plus_bjet1 = baseLepton->mom() + mostBjetLike[1]->mom();

      double pa_a[3] = {0, lepton_plus_bjet0.px(), lepton_plus_bjet0.py()};
      double pb_a[3] = {80, mostBjetLike[1]->mom().px(), mostBjetLike[1]->mom().py()};
      double pmiss_a[3] = {0, metVec.px(), metVec.py()};
      double mn_a = 0.;

      mt2_bisect::mt2 mt2_event_a;

      mt2_event_a.set_momenta(pa_a, pb_a, pmiss_a);
      mt2_event_a.set_mn(mn_a);

      double mt2a = mt2_event_a.get_mt2();

      double pa_b[3] = {0, lepton_plus_bjet1.px(), lepton_plus_bjet1.py()};
      double pb_b[3] = {80, mostBjetLike[0]->mom().px(), mostBjetLike[0]->mom().py()};
      double pmiss_b[3] = {0, metVec.px(), metVec.py()};
      double mn_b = 0.;

      mt2_bisect::mt2 mt2_event_b;

      mt2_event_b.set_momenta(pa_b, pb_b, pmiss_b);
      mt2_event_b.set_mn(mn_b);
      double mt2b = mt2_event_b.get_mt2();
      return min(mt2a, mt2b);
    }

    double Analysis_ATLAS_13TeV_1LEPStop_36invfb::calcMetPerp(const HEPUtils::P4& topChi2, const HEPUtils::P4& top1)
    {
      HEPUtils::P4 ttbarP4 = topChi2 + top1;

      TLorentzVector ttbar, top1Rest, metRest;
      ttbar.SetPxPyPzE(ttbarP4.px(), ttbarP4.py(), ttbarP4.pz(), ttbarP4.E());
      top1Rest.SetPxPyPzE(top1.px(), top1.py(), top1.pz(), top1.E());
      metRest.SetPxPyPzE(metVec.px(), metVec.py(), metVec.pz(), metVec.E());

      TVector3 boostVector = TVector3(-ttbar.Px() / ttbar.E(), -ttbar.Py() / ttbar.E(), -ttbar.Pz() / ttbar.E());
      top1Rest.Boost(boostVector);
      metRest.Boost(boostVector);
      return metRest.Vect().XYvector().Norm(top1Rest.Vect().XYvector()).Mod();
    }

    double Analysis_ATLAS_13TeV_1LEPStop_36invfb::calcMT2Tau(vector<HEPUtils::Particle*> baselineTaus, HEPUtils::Particle* baseLepton)
    {
      double mT2Tau = 120.0;
      if (!baselineTaus.empty()) {
        double pa_tau[3] = {0, baselineTaus[0]->mom().px(), baselineTaus[0]->mom().py()};
        double pb_tau[3] = {0, baseLepton->mom().px(), baseLepton->mom().py()};
        double pmiss_tau[3] = {0, metVec.px(), metVec.py()};
        double mn_tau = 0.;
        mt2_bisect::mt2 mt2_event_tau;
        mt2_event_tau.set_momenta(pa_tau, pb_tau, pmiss_tau);
        mt2_event_tau.set_mn(mn_tau);
        mT2Tau = mt2_event_tau.get_mt2();
      }
      return mT2Tau;
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::assignAbsDPhiJMet(double absDPhiJMet[4])
    {
      size_t i;
      for (i = 0; i < 4; i++)
      {
        i < nJets ? absDPhiJMet[i] = fabs(signalJets[i]->mom().deltaPhi(metVec)) : NAN;
      }
    }


    inline void Analysis_ATLAS_13TeV_1LEPStop_36invfb::addBDTVariable(
      TMVA::Reader* reader,
      const string& varName,
      double value)
    {
      reader->AddVariable(varName, new float(value));
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::analyze(const HEPUtils::Event* event)
    {
      // must call initVars here to reset some crucial variables that are updated on each analyze call
      initVars();

      // Missing energy
      metVec = event->missingmom();
//      Met = event->met();
      metVec.setPE(-metVec.px(), -metVec.py(), -metVec.pz(), metVec.E());
      Met = metVec.pT();

      // Baseline lepton objects
      vector<HEPUtils::Particle*> baselineElectrons = filterPtEta(event->electrons(), 5., 2.47);
      vector<HEPUtils::Particle*> baselineMuons = filterPtEta(event->muons(), 4., 2.7);
      vector<HEPUtils::Particle*> baselineTaus = filterPtEta(event->taus(), 20.0, 2.5);
      ATLAS::applyTauEfficiencyR1(baselineTaus);

      // Get b jets
      vector<HEPUtils::Jet*> bJets, nonBJets;
      getBJets(event, &bJets, &nonBJets);

      // Overlap removal
      vector<HEPUtils::Particle*> signalElectrons, signalSoftElectrons, signalMuons, signalSoftMuons;
      vector<HEPUtils::Particle*> signalSoftLeptons, electronsForVeto, muonsForVeto;

      // Note: use paper description instead of code snippet
      // This is not identical to the overlap removal in the paper
      // Probably good enough though
      LeptonLeptonOverlapRemoval(&baselineMuons, &baselineElectrons, 0.01); // mimics shared track requirement
      JetLeptonOverlapRemoval(&nonBJets, &baselineElectrons, 0.2);
      LeptonJetOverlapRemoval(&baselineElectrons, &nonBJets);
      LeptonJetOverlapRemoval(&baselineElectrons, &bJets);
      LeptonJetOverlapRemoval(&baselineMuons, &nonBJets);
      LeptonJetOverlapRemoval(&baselineMuons, &bJets);
      LeptonLeptonOverlapRemoval(&baselineTaus, &baselineElectrons, 0.1);

      // Put the sorted baseline electrons and muons into the baseline leptons vector
      baselineLeptons.insert(baselineLeptons.end(), baselineElectrons.begin(), baselineElectrons.end());
      baselineLeptons.insert(baselineLeptons.end(), baselineMuons.begin(), baselineMuons.end());

      vector<HEPUtils::Jet*> signalBJets = getSignalTypeJets(bJets, 25., 2.5, true);
      vector<HEPUtils::Jet*> signalNonBJets = getSignalTypeJets(nonBJets, 25., 2.5, false);

      /* ensure object collections to be pT sorted */
      sort(signalJets.begin(), signalJets.end(), sortByPT_1l);

      // Note that the isolation requirements and tight selection are currently missing
      fillSignalParticles(&baselineElectrons, &signalSoftElectrons, &signalElectrons, &signalSoftLeptons, 25.0);
      fillSignalParticles(&baselineMuons, &signalSoftMuons, &signalMuons, &signalSoftLeptons, 25.0);

      // We now have the signal electrons, muons, jets and b jets- move on to the analysis

      nJets = signalJets.size();
      nBJets = signalBJets.size();

      // Minimal event selection
      // if failed then we should add to our cut flows and bail out of the analysis for this event
      if (!(Met > 100. &&
            baselineLeptons.size() == 1 &&
            (signalSoftLeptons.size() == 1 || signalLeptons.size() == 1) &&
            nJets > 1))
      {
        return;
      }

      HEPUtils::Particle* baseLepton = baselineLeptons[0];

      vector<HEPUtils::Jet*> mostBjetLike, signalNotBjetLike, signalNotBjet;

      vector<int> bJetNums = bJetSelection();
      mostBjetLike.push_back(signalJets[bJetNums[0]]);
      mostBjetLike.push_back(signalJets[bJetNums[1]]);

      getSignalNotBjets(&signalNotBjet, &signalNotBjetLike, bJetNums);

      /* ensure object collections to be pT sorted */
      sort(signalJets.begin(), signalJets.end(), sortByPT_1l);

      if (!baselineTaus.empty())
      {
        sort(baselineTaus.begin(), baselineTaus.end(), sortByPT_1l_Particle);
      }

      // Now make a collection to hold the JER for each jet
      // Have obtained the values from Matthias' BuckFast code
      // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2015-017/
      // Parameterisation can be still improved, but eta dependence is minimal
      const vector<double> binedges_eta = {0, 10.};
      const vector<double> binedges_pt = {0, 50., 70., 100., 150., 200., 1000., 10000.};
      const vector<double> JetsJER = {0.145, 0.115, 0.095, 0.075, 0.07, 0.05, 0.04};
      static HEPUtils::BinnedFn2D<double> _resJets2D(binedges_eta, binedges_pt, JetsJER);

      for (HEPUtils::Jet* jet : signalJets)
      {
        signalJER.push_back(_resJets2D.get_at(jet->abseta(), jet->pT()));
      }

      HtSigMiss = calcHtSigMiss();

      double absDPhiJMet[4];
      assignAbsDPhiJMet(absDPhiJMet);

      absDPhiJMet0 = absDPhiJMet[0];
      absDPhiJMet1 = absDPhiJMet[1];

      double absDPhiJiMet = getAbsDPhiJiMet(absDPhiJMet);

      mT = calcMT_1l(baseLepton->mom(), metVec);

      dPhiMetLep = fabs(metVec.deltaPhi(baseLepton->mom()));

      // Calculate MT2 tau using the leading tau in the event
      mT2Tau = calcMT2Tau(baselineTaus, baseLepton);

      pTLepOverMet = baseLepton->pT() / Met;
      bool preselHighMet = Met > 230 && mT > 30;
      bool preselLowMet =
        baseLepton->pT() > 27 &&
        !signalBJets.empty() &&
        signalJets[0]->pT() > 50. &&
        Met > 100 &&
        mT > 90;

      // Apply tight selection if lepton is an electron
      // Am using same selection as 8 TeV (probably needs updating)
      // Note that we have already applied a 1 lepton cut
      if (baselineElectrons.size() == 1 && baselineMuons.empty()) {
        vector<HEPUtils::Particle*> tightElectrons;
        tightElectrons.push_back(baselineElectrons[0]);
        ATLAS::applyTightIDElectronSelection(tightElectrons);
        preselLowMet = preselLowMet && (tightElectrons.size() == 1);
      }

      // Now calculate amT2 using two different assignments of the b jets and the leptons
      amT2 = calcAmt2(baseLepton, mostBjetLike);
      dRbl = baseLepton->mom().deltaR_eta(mostBjetLike[0]->mom());

      /* Reconstruct top by a chi2 based method */
      vector<int> jetComb = getJetComb(bJetNums, mostBjetLike, _resJets2D);
      HEPUtils::P4 topChi2 = signalJets[jetComb[0]]->mom() + signalJets[jetComb[1]]->mom() + mostBjetLike[jetComb[2]]->mom();
      HEPUtils::P4 top1 = baseLepton->mom() + (jetComb[2] == 0 ? mostBjetLike[1]->mom() : mostBjetLike[0]->mom());

      // calculate MetPerp
      MetPerp = calcMetPerp(topChi2, top1);

      // Do RestFrames Calculations
      doRJRCalculations();

      // Now we have to do the fancy jet reclustering to get reconstructed W and top particles
      // signalNotBjet + mostBjetLike is inconsistent but bjets are not used anyway
      HEPUtils::P4 WRecl = reclusteredParticle(signalNotBjet, mostBjetLike, mW, false);
      HEPUtils::P4 topRecl = reclusteredParticle(signalNotBjetLike, mostBjetLike, 175., true);

      if (nBJets > 0 && nJets > 3 && preselHighMet)
      {
        topReclM = topRecl.m();
      }

      WReclM = WRecl.m();

      // Should now be ready to do signal selections

      // need to assign these now because we have to check that the signal jets are actually there at each index
      assignJetValues(signalBJets);

      // use this vector for checking multiple jet pTs at once
      vector<double> signalJetPts = {jetPt1, jetPt2, jetPt3, jetPt4, jetPt5};

      double dPhiTTBar = dPhi(top1.phi(), topChi2.phi());

      if (preselLowMet)
      {
        varResults["tN_diag_med_preMTopChi2"].push_back(topChi2.m());
        varResults["tN_diag_med_MET"].push_back(Met);
        varResults["tN_diag_med_mT"].push_back(mT);
        varResults["tN_diag_med_htSigMiss"].push_back(HtSigMiss);
        varResults["tN_diag_med_dRbl"].push_back(dRbl);
        varResults["tN_diag_med_dPhiTTbar"].push_back(dPhiTTBar);
        varResults["dPhiHadTopMet"].push_back(dPhi(topChi2.phi(), metVec.phi()));
        varResults["tN_diag_med_nJets"].push_back(nJets);
        varResults["jet_pt[2]"].push_back(signalJets[2]->pT());
        varResults["jet_pt[3]"].push_back(signalJets[3]->pT());
      }

      double gev = 1000;

      // non-soft lepton selections
      // All of the numbers here are specified in the tables 6 - 11 in the paper. An overview is given in table 5
      if (signalLeptons.size() == 1)
      {
        // tN_med - table 6
        acceptSignalRegion(
          tN_med,
          nJets > 3 &&
          nBJets > 0 &&
          preselLowMet &&
          jetPt1 > 60 &&
          jetPt2 > 50 &&
          jetPt4 > 40 &&
          Met > 250 &&
          MetPerp > 230 &&
          HtSigMiss > 14 &&
          mT > 160 &&
          amT2 > 175 &&
          topReclM > 150 &&
          dRbl < 2.0 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          mT2Tau > 80);

        // tN_high - table 6
        acceptSignalRegion(
          tN_high,
          nJets > 3 &&
          nBJets > 0 &&
          preselHighMet &&
          // the array of numbers here is from the jet pT requirements from the table
          checkJetPt({100, 80, 50, 30}, signalJetPts) &&
          Met > 550 &&
          HtSigMiss > 27 &&
          mT > 160 &&
          amT2 > 175 &&
          topReclM > 130 &&
          dRbl < 2.0 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          mT2Tau > 80);

        // bWN - table 8
        acceptSignalRegion(
          bWN,
          nJets > 3 &&
          nBJets > 0 &&
          preselHighMet &&
          jetPt1 > 50 &&
          Met > 300 &&
          mT > 130 &&
          amT2 < 110 &&
          dPhiMetLep < 2.5 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          mT2Tau > 80);

        // bC2x_diag - table 9
        acceptSignalRegion(
          bC2x_diag,
          nJets > 3 &&
          nBJets > 1 &&
          preselHighMet &&
          jetPt3 > 75 &&
          jetPt4 > 30 &&
          jetPt2 > 30 &&
          Met > 230 &&
          HtSigMiss > 13 &&
          mT > 180 &&
          amT2 > 175 &&
          absDPhiJMet0 > 0.7 &&
          absDPhiJMet1 > 0.7 &&
          WReclM > 50 &&
          mT2Tau > 80);

        // bC2x_med - table 9
        acceptSignalRegion(
          bC2x_med,
          nJets > 3 &&
          nBJets > 1 &&
          preselHighMet &&
          checkJetPt({200, 140}, signalJetPts) &&
          bJetPt2 > 140 &&
          Met > 230 &&
          HtSigMiss > 10 &&
          mT > 120 &&
          amT2 > 300 &&
          absDPhiJMet0 > 0.9 &&
          absDPhiJMet1 > 0.9 &&
          WReclM > 50 &&
          mT2Tau > 80);

        // bCbv - table 9
        acceptSignalRegion(
          bCbv,
          nJets > 1 &&
          nBJets == 0 &&
          preselHighMet &&
          checkJetPt({120, 80}, signalJetPts) &&
          Met > 360 &&
          HtSigMiss > 16 &&
          mT > 200 &&
          absDPhiJMet0 > 2.0 &&
          absDPhiJMet1 > 0.8 &&
          WReclM >= 70 &&
          WReclM <= 100 &&
          dPhiMetLep > 1.2 &&
          baseLepton->pT() > 60);

        // DM_low - table 11
        acceptSignalRegion(
          DM_low,
          nJets > 3 &&
          nBJets > 0 &&
          preselHighMet &&
          checkJetPt({120, 85, 65}, signalJetPts) &&
          bJetPt1 > 60 &&
          Met > 320 &&
          mT > 170 &&
          HtSigMiss > 14 &&
          amT2 > 160 &&
          topReclM > 130 &&
          dPhiMetLep > 1.2 &&
          absDPhiJiMet > 1.0 &&
          mT2Tau > 80);

        // DM_high - table 11
        acceptSignalRegion(
          DM_high,
          nJets > 3 &&
          nBJets > 0 &&
          preselHighMet &&
          checkJetPt({125, 75, 65}, signalJetPts) &&
          Met > 380 &&
          mT > 225 &&
          amT2 > 190 &&
          topReclM > 130 &&
          dPhiMetLep > 1.2 &&
          absDPhiJiMet > 1.0);


        // tN_diag_low - table 7
        if (preselLowMet &&
            nJets >= 4 &&
            nBJets >= 1 &&
            jetPt1 > 120 &&
            Met > 100 &&
            mT > 90 &&
            absDPhiJMet0 > 0.4 &&
            absDPhiJMet1 > 0.4)
        {
          HEPUtils::P4 topHadronic, bLepton, nu, ttbar, rjet;
          chi2Selection(&topHadronic, &bLepton, mostBjetLike, signalNotBjetLike);

          HEPUtils::P4 signalLepton = *signalLeptons[0];

          HEPUtils::P4 particleSum = topHadronic + bLepton + signalLepton;

          const double alpha = 27.0 / 200.0;
          auto nuPtComponent = [&alpha](double ptMiss, double ptSum) { return (1.0 - alpha) * ptMiss - alpha * ptSum; };

          double nuPx = nuPtComponent(metVec.px(), particleSum.px());
          double nuPy = nuPtComponent(metVec.py(), particleSum.py());
          double nuPhi = atan2(nuPy, nuPx);
          double nuPt = sqrt(sqr(nuPx) + sqr(nuPy));
          // in the internal document ATL-COM-PHYS-2017-214 it states that we should use the signal lepton's
          // pseudorapidity (eta) for the neutrino so that we can solve for the leptonic top mass
          nu.setEtaPhiMPt(signalLepton.eta(), nuPhi, 0.0, nuPt);

          HEPUtils::P4 topLeptonic = bLepton + signalLepton + nu;

          ttbar = particleSum + metVec; // under SM hypothesis

          rjet.setPE(-1.0 * ttbar.px(), -1.0 * ttbar.py(), -1.0 * ttbar.pz(), ttbar.E());

          // calculate the rest of the variables we need for the BDT
          double dPhiLeptonNu = signalLepton.deltaPhi(nu);
          double mTalpha = sqrt((1.0 - cos(dPhiLeptonNu)) * 2.0 * signalLepton.pT() * nu.pT());
          double DMTalpha = mT - mTalpha;

          double MTopLep200 = topLeptonic.m();
          double MTopHadronic = topHadronic.m();
          double DPhiLepNu200 = fabs(dPhiLeptonNu);
          double DPhiRJetLep  = fabs(signalLepton.deltaPhi(rjet));

          varResults["preBLepMass"].push_back(bLepton.m());
          varResults["preTopHadronic"].push_back(MTopHadronic);
          varResults["preRjetPt"].push_back(rjet.pT());
          varResults["preDPhiRjetLep"].push_back(DPhiRJetLep);
          varResults["preMTopLep200"].push_back(MTopLep200);
          varResults["preDMTAlpha"].push_back(DMTalpha);
          varResults["preDPhiLepNu200"].push_back(DPhiLepNu200);
          varResults["preMet"].push_back(Met);
          varResults["preMt"].push_back(mT);
          varResults["preNjets"].push_back(nJets);
          varResults["preNBjets"].push_back(nBJets);

          if (rjet.pT() > 400 && DPhiRJetLep > 1.0)
          {
            varResults["DMTAlpha"].push_back(DMTalpha);
            varResults["tN_diag_low_MET"].push_back(Met);
            varResults["tN_diag_low_mT"].push_back(mT);
            varResults["topHadronicMass"].push_back(MTopHadronic);
            varResults["topLeptonic200Mass"].push_back(MTopLep200);
            varResults["dPhiLepNu200"].push_back(DPhiLepNu200);
            varResults["dPhiRjetLep"].push_back(DPhiRJetLep);

            string methodName = "ATLAS_13TeV_1LEPStop_36invfbtN_diag_low_BDT";
            TMVA::Reader* reader = new TMVA::Reader("Silent");
            addBDTVariable(reader, "MET", Met * gev);
            addBDTVariable(reader, "MT", mT * gev);
            addBDTVariable(reader, "dMT200", DMTalpha * gev);
            addBDTVariable(reader, "m_tophad", MTopHadronic * gev);
            addBDTVariable(reader, "m_toplep200", MTopLep200 * gev);
            addBDTVariable(reader, "dphi_lep_nu200", DPhiLepNu200);
            addBDTVariable(reader, "dphi_rjet_lep", DPhiRJetLep);
            reader->BookMVA(methodName, BDTWeightsLocation + "BDT_low" + BDTWeightsFileNum + ".weights.xml");
            double tN_diag_low_BDTScore = reader->EvaluateMVA(methodName);
            BDTLowScores.push_back(tN_diag_low_BDTScore);
            acceptSignalRegion(tN_diag_low, tN_diag_low_BDTScore > 0.55);
          }
        }

        // tN_diag_med - table 7
        if (preselLowMet &&
            nJets >= 4 &&
            nBJets >= 1 &&
            jetPt1 > 100 &&
            jetPt2 > 50 &&
            Met > 120 &&
            mT > 120 &&
            absDPhiJMet0 > 0.4 &&
            absDPhiJMet1 > 0.4 &&
            mT2Tau > 80)
        {
          varResults["mTopChi2"].push_back(topChi2.m());

          string methodName = "ATLAS_13TeV_1LEPStop_36invfbtN_diag_med_BDT";
          TMVA::Reader* reader = new TMVA::Reader("Silent");
          addBDTVariable(reader, "met", Met * gev);
          addBDTVariable(reader, "mt", mT * gev);
          addBDTVariable(reader, "m_top_chi2", topChi2.m());
          addBDTVariable(reader, "ht_sig", HtSigMiss);
          addBDTVariable(reader, "dr_bjet_lep", dRbl);
          addBDTVariable(reader, "dphi_ttbar", dPhiTTBar);
          addBDTVariable(reader, "dphi_hadtop_met", dPhi(topChi2.phi(), metVec.phi()));
          addBDTVariable(reader, "n_jet", nJets);
          addBDTVariable(reader, "jet_pt[2]", signalJets[2]->pT() * gev);
          addBDTVariable(reader, "jet_pt[3]", signalJets[3]->pT() * gev);
          reader->AddSpectator("sf_total", new float(1.0));
          reader->AddSpectator("xs_weight", new float(1.0));
          reader->AddSpectator("weight", new float(1.0));
          reader->BookMVA(methodName, BDTWeightsLocation + "BDT_med" + BDTWeightsFileNum + ".weights.xml");
          double tN_diag_med_BDTScore = reader->EvaluateMVA(methodName);
          BDTMedScores.push_back(tN_diag_med_BDTScore);
          acceptSignalRegion(tN_diag_med, tN_diag_med_BDTScore > 0.75);
        }

        // tN_diag_high - table 7
        if (preselHighMet &&
            nJets >= 5 &&
            nBJets >= 1 &&
            Met > 230 &&
            mT > 120 &&
            rISR > 0.4)
        {
          if (checkJetPt({25, 25, 25, 25}, {jetPt2, jetPt3, jetPt4, jetPt5}))
          {
            varResults["mTS"].push_back(mTS);
            varResults["rISR"].push_back(rISR);
            varResults["dPhiISR"].push_back(dphiISRI);
            varResults["tN_diag_high_mT"].push_back(mT);
            varResults["tN_diag_high_dPhiTTbar"].push_back(dPhiTTBar);
            varResults["tN_diag_high_dRbl"].push_back(dRbl);
            varResults["jetPt3"].push_back(jetPt3);
            varResults["jetPt4"].push_back(jetPt4);
            varResults["tN_diag_high_mTopChi2"].push_back(topChi2.m());
            varResults["njv"].push_back(njv);
            varResults["tN_diag_high_nJets"].push_back(nJets);
            varResults["tN_diag_high_MET"].push_back(Met);
          }
          string methodName = "ATLAS_13TeV_1LEPStop_36invfbtN_diag_high_BDT";
          TMVA::Reader* reader = new TMVA::Reader("Silent");
          addBDTVariable(reader, "mt", mT * gev);
          addBDTVariable(reader, "dphi_ttbar", dPhiTTBar);
          addBDTVariable(reader, "dr_bjet_lep", dRbl);
          addBDTVariable(reader, "jetPt3", jetPt3 * gev);
          addBDTVariable(reader, "jetPt4", jetPt4 * gev);
          addBDTVariable(reader, "m_top_chi2", topChi2.m());
          addBDTVariable(reader, "RJR_RISR", rISR);
          addBDTVariable(reader, "RJR_MS", mTS);
          addBDTVariable(reader, "RJR_dphiISRI", dphiISRI);
          addBDTVariable(reader, "RJR_NjV", njv);
          reader->BookMVA(methodName, BDTWeightsLocation + "BDT_high" + BDTWeightsFileNum + ".weights.xml");
          double tN_diag_high_BDTScore = reader->EvaluateMVA(methodName);
          BDTHighScores.push_back(tN_diag_high_BDTScore);
          acceptSignalRegion(tN_diag_high, tN_diag_high_BDTScore > 0.8);
        }
      }

      // Soft-lepton selections
      for (HEPUtils::Jet* jet : signalBJets)
      {
        double dPhi_tmp = fabs(jet->mom().deltaPhi(metVec));
        if (dPhi_tmp < minDPhiMetBJet)
        {
          minDPhiMetBJet = dPhi_tmp;
        }
      }

      if (nBJets > 1)
      {
        dRbb = mostBjetLike[0]->mom().deltaR_eta(mostBjetLike[1]->mom());
      }


      Wpt = NAN;
      if (signalSoftLeptons.size() == 1)
      {
        preselSoftLep = Met > 230;

        // Apply tight selection if lepton is an electron
        // Am using same selection as 8 TeV (probably needs updating)
        // Note that we have already applied a 1 lepton cut
        if (signalSoftElectrons.size() == 1 && signalSoftMuons.empty())
        {
          vector<HEPUtils::Particle*> tightElectrons;
          tightElectrons.push_back(signalSoftElectrons[0]);
          ATLAS::applyTightIDElectronSelection(tightElectrons);
          preselSoftLep = preselSoftLep && (tightElectrons.size() == 1);
        }

        Wpt = (signalSoftLeptons[0]->mom() + metVec).pT();

        // bffN - table 8
        acceptSignalRegion(
          bffN,
          nJets > 1 &&
          nBJets > 0 &&
          preselSoftLep &&
          jetPt1 > 400 &&
          Met > 300 &&
          mT < 160 &&
          pTLepOverMet < 0.02 &&
          minDPhiMetBJet < 1.5 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          topReclM < 150 &&
          !jet1Btag);

        // bCsoft_diag - table 10
        acceptSignalRegion(
          bCsoft_diag,
          nJets > 1 &&
          nBJets > 0 &&
          preselSoftLep &&
          jetPt1 > 400 &&
          Met > 300 &&
          mT < 50 &&
          pTLepOverMet < 0.02 &&
          minDPhiMetBJet < 1.5 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          topReclM < 150 &&
          !jet1Btag);

        // bCsoft_med - table 10
        acceptSignalRegion(
          bCsoft_med,
          nJets > 2 &&
          nBJets > 1 &&
          preselSoftLep &&
          checkJetPt({120, 60, 40}, signalJetPts) &&
          checkJetPt({120, 60}, signalJetPts) &&
          Met > 230 &&
          mT < 160 &&
          pTLepOverMet < 0.03 &&
          amT2 > 200 &&
          minDPhiMetBJet > 0.8 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          Wpt > 400);

        // bCsoft_high - table 10
        acceptSignalRegion(
          bCsoft_high,
          nJets > 1 &&
          nBJets > 1 &&
          preselSoftLep &&
          jetPt2 > 100 &&
          bJetPt2 > 100 &&
          Met > 230 &&
          mT < 160 &&
          pTLepOverMet < 0.03 &&
          amT2 > 300 &&
          minDPhiMetBJet > 0.4 &&
          absDPhiJMet0 > 0.4 &&
          absDPhiJMet1 > 0.4 &&
          Wpt > 500 &&
          dRbb > 0.8);
      }

    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::acceptSignalRegion(SignalRegion signalRegion, bool isAccepted)
    {
      if (isAccepted)
      {
        signalRegionCounts[SignalRegionStrings[signalRegion]]++;
      }
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::saveCutFlow()
    {
      double scale_by = 1.0;
      cout << "SAVE_START:cuts" << endl;
      cout << "CUT;RAW;SCALED;%" << endl;
      double initialCut = cutFlowVector[0];
      double thisCut;
      for (unsigned j = 0; j < NCUTS; j++) {
        thisCut = cutFlowVector[j];
        cout << cutFlowVector_str[j].c_str() << ";"
             << thisCut << ";"
             << thisCut * scale_by << ";"
             << 100. * thisCut / initialCut
             << endl;
      }
      cout << "SAVE_END" << endl;
    }



    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::save1DVector(vector<double> vec, const string &fileName, int numEvents)
    {
      std::ofstream outputFile;
      std::string N = std::to_string(numEvents);
      std::string outFileName = "/home/dsteiner/madgraph/analysis/data/" + fileName + "_" + N + ".csv";
      cout << outFileName << endl;
      outputFile.open(outFileName);
      for (double value : vec) {
        outputFile << value << endl;
      }
      outputFile.close();
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::collect_results(int numEvents)
    {
      // For debugging purposes
//      saveCutFlow();
      save1DVector(BDTLowScores, "tN_diag_low", numEvents);
      save1DVector(BDTMedScores, "tN_diag_med", numEvents);
      save1DVector(BDTHighScores, "tN_diag_high", numEvents);
      for (pair<string, vector<double>> entry : varResults)
      {
        save1DVector(entry.second, entry.first, numEvents);
      }
    }

    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::initVars()
    {
      baselineLeptons = {};
      signalLeptons = {};
      signalJets = {};
      signalJER = {};
      nJets = 0;
      nBJets = 0;
      Met = 0.0;
      absDPhiJMet0 = 0.0;
      absDPhiJMet1 = 0.0;
      mT2Tau = 0.0;
      jetPt1 = 0.0;
      jetPt2 = 0.0;
      jetPt3 = 0.0;
      jetPt4 = 0.0;
      jetPt5 = 0.0;
      MetPerp = 0.0;
      HtSigMiss = 0.0;
      mT = 0.0;
      amT2 = 0.0;
      dRbl = 0.0;
      topReclM = 0.0;
      dPhiMetLep = 0.0;
      minDPhiMetBJet = DBL_MAX;
      pTLepOverMet = 0.0;
      bJetPt1 = 0.0;
      bJetPt2 = 0.0;
      WReclM = 0.0;
      Wpt = 0.0;
      dRbb = 0.0;
      mW = 80.0;
      mTop = 170.0;
      njv = 0;
      rISR = 0.0;
      dphiISRI = 0.0;
      mTS = 0.0;
      preselSoftLep = false;
      // initialize this to true because we usually check if its not true
      jet1Btag = true;
    }


    void Analysis_ATLAS_13TeV_1LEPStop_36invfb::clear()
    {
      fill(cutFlowVector.begin(), cutFlowVector.end(), 0);
      BDTLowScores.clear();
      BDTMedScores.clear();
      BDTHighScores.clear();
      for (pair<string, int> entry : signalRegionCounts)
      {
        signalRegionCounts[entry.first] = 0;
      }
      for (pair<string, vector<double>> entry : varResults)
      {
        varResults[entry.first] = vector<double>({});
      }
    }
  }
}

void tev13_analysis(const char *inputFile, const char *outputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  double massElectron = 0.000510998902;
  double massMuon = 0.105658389;

  Gambit::ColliderBit::Analysis_ATLAS_13TeV_1LEPStop_36invfb analysis = Gambit::ColliderBit::Analysis_ATLAS_13TeV_1LEPStop_36invfb();
  analysis.clear();
  //------------------------------------------------------------------------------
  // Loop over all events and write to output file
  int numEvents = 0;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    if (entry == 10000)
    {
      exit(0);
    }
    numEvents++;
    if (std::fmod(entry + 1, 10000) == 0)
    {
      cout << "Entry number: " << entry + 1 << endl;
    }
    treeReader->ReadEntry(entry);
    //------------------------------------------------------------------------------
    //Convert to the form used in GAMBIT code
    HEPUtils::Event eventPRE;
    Candidate *candidate;
    eventPRE.clear();

    HEPUtils::Particle *recoParticle;
    HEPUtils::Jet *recoJet;
    // Delphes particle arrays: Post-Detector Sim
    //    MISSING ET:
    MissingET *met = (MissingET *) branchMissingET->At(0);
    TLorentzVector momentum;// = candidate->Momentum;
    momentum.SetPtEtaPhiM(met->MET, met->Eta, met->Phi, 0);
    eventPRE.set_missingmom(HEPUtils::P4::mkXYZM(-1 * momentum.Px(), -1 * momentum.Py(), 0., 0.));

    // Delphes particle arrays: Post-Detector Sim
    //    PHOTONS:
    for (int i = 0; i < branchPhoton->GetEntriesFast(); i++) {
      Photon *candidate = (Photon *) branchPhoton->At(i);
      TLorentzVector momentum;// = candidate->Momentum;
      momentum.SetPtEtaPhiM(candidate->PT, candidate->Eta, candidate->Phi, 0.);
      recoParticle = new HEPUtils::Particle(HEPUtils::P4::mkXYZM(momentum.Px(), momentum.Py(), momentum.Pz(), 0.),
                                            MCUtils::PID::PHOTON);
      recoParticle->set_prompt(true);
      eventPRE.add_particle(recoParticle);
    }

    // Delphes particle arrays: Post-Detector Sim
    //    ELECTRONS:
    for (int i = 0; i < branchElectron->GetEntriesFast(); i++) {
      Electron *candidate = (Electron *) branchElectron->At(i);
      TLorentzVector momentum;// = candidate->Momentum;
      momentum.SetPtEtaPhiM(candidate->PT, candidate->Eta, candidate->Phi, massElectron);
      recoParticle = new HEPUtils::Particle(
        HEPUtils::P4::mkXYZM(momentum.Px(), momentum.Py(), momentum.Pz(), 0.000510998902),
        -HEPUtils::sign(candidate->Charge) * MCUtils::PID::ELECTRON);
      recoParticle->set_prompt(true);
      eventPRE.add_particle(recoParticle);
    }

    // Delphes particle arrays: Post-Detector Sim
    //    MUONS:
    for (int i = 0; i < branchMuon->GetEntriesFast(); i++) {
      Muon *candidate = (Muon *) branchMuon->At(i);
      TLorentzVector momentum;// = candidate->Momentum;
      momentum.SetPtEtaPhiM(candidate->PT, candidate->Eta, candidate->Phi, massMuon);
      recoParticle = new HEPUtils::Particle(
        HEPUtils::P4::mkXYZM(momentum.Px(), momentum.Py(), momentum.Pz(), 0.105658389),
        -HEPUtils::sign(candidate->Charge) * MCUtils::PID::MUON);
      recoParticle->set_prompt(true);
      eventPRE.add_particle(recoParticle);
    }

    // Delphes particle arrays: Post-Detector Sim
    //    JETS and TAUS:
    for (int i = 0; i < branchJet->GetEntriesFast(); i++) {
      Jet *candidate = (Jet *) branchJet->At(i);
      TLorentzVector momentum;// = candidate->Momentum;
      momentum.SetPtEtaPhiM(candidate->PT, candidate->Eta, candidate->Phi, candidate->Mass);
      if (candidate->TauTag) {
        recoParticle = new HEPUtils::Particle(HEPUtils::P4::mkXYZM(momentum.Px(), momentum.Py(), momentum.Pz(), 1e-6),
                                              -HEPUtils::sign(candidate->Charge) * MCUtils::PID::TAU);
        recoParticle->set_prompt(true);
        eventPRE.add_particle(recoParticle);
        //continue;
      } else {
        /// @todo Should the jet mass be assigned properly rather than set as microscopic?
        recoJet = new HEPUtils::Jet(HEPUtils::P4::mkXYZM(momentum.Px(), momentum.Py(), momentum.Pz(), 1e-6),
                                    candidate->BTag);
        eventPRE.add_jet(recoJet);
      }
    }

    const HEPUtils::Event* event = &eventPRE;

    asymm_mt2_lester_bisect::disableCopyrightMessage();

    analysis.analyze(event);
  }
  analysis.collect_results(numEvents);
}
