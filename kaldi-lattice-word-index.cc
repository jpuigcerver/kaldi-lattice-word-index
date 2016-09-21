#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

namespace fst {

void ConvertLatticeWeight(const kaldi::CompactLatticeWeight& iw,
                          StdArc::Weight* ow) {
  KALDI_ASSERT(ow != NULL);
  *ow = StdArc::Weight(iw.Weight().Value1() + iw.Weight().Value2());
}

void ConvertLatticeWeight(const kaldi::CompactLatticeWeight& iw,
                          LogArc::Weight* ow) {
  KALDI_ASSERT(ow != NULL);
  *ow = LogArc::Weight(iw.Weight().Value1() + iw.Weight().Value2());
}

template <typename F1, typename F2>
void ConvertFst(const F1& f1, F2* f2) {
  typedef typename F1::Arc A1;
  typedef typename F2::Arc A2;
  fst::ArcMap(f1, f2, fst::WeightConvertMapper<A1, A2>());
}

}  // namespace fst

namespace kaldi {

void AddInsPenToLattice(BaseFloat penalty, CompactLattice *lat) {
  typedef typename CompactLattice::Arc Arc;
  typedef typename CompactLattice::Weight::W Weight;
  for (int32 state = 0; state < lat->NumStates(); ++state) {
    for (fst::MutableArcIterator<CompactLattice> aiter(lat, state);
         !aiter.Done(); aiter.Next()) {
      Arc arc(aiter.Value());
      if (arc.olabel != 0) {
        Weight weight = arc.weight.Weight();
        weight.SetValue1(weight.Value1() + penalty);
        arc.weight.SetWeight(weight);
        aiter.SetValue(arc);
      }
    }
  }
}

}  // namespace kaldi

// Given a fst that represents sequences of characters in its paths, creates
// a fst that accepts words (determined by sequences of characters in between
// special labels, separators) which where part of the original fst.
template <typename Arc>
void WordsFst(fst::MutableFst<Arc>* fst,
              const std::unordered_set<typename Arc::Label>& separators) {
  typedef fst::MutableFst<Arc> Fst;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  KALDI_ASSERT(fst != NULL);
  if (fst->Start() == fst::kNoStateId) return;

  // Compute forward and backward scores for each state
  std::vector<Weight> fw, bw;
  fst::ShortestDistance<Arc>(*fst, &fw, false);
  fst::ShortestDistance<Arc>(*fst, &bw, true);
  const float total_cost = bw[fst->Start()].Value();

  // New final state
  const StateId sFinal = fst->AddState();

  // Convert fst to accept words in the original fst
  std::vector<Arc> new_arcs_from_init;
  for (fst::StateIterator<Fst> siter(*fst); !siter.Done(); siter.Next()) {
    const StateId s = siter.Value();
    if (s == sFinal) continue;
    std::vector<Arc> new_arcs;  // New arcs from state s
    // Add arc to the new (unique) final state, and make s not final
    if (fst->Final(s) != Weight::Zero()) {
      new_arcs.push_back(Arc(0, 0, fst->Final(s), sFinal));
      fst->SetFinal(s, Weight::Zero());
    }
    // Traverse current arcs and remove arcs with separator labels.
    // For each remove arc two epsilon arcs are added:
    //   - one from the current state (s) to the final state with cost =
    //     arc.weight * bw[arc.nextstate]. This is because the node s is
    //     the final node for a word.
    //   - one from the initial state to arc.nextstate, with cost =
    //     arc.weight * forward[s]. This is because arc.nextstate is the
    //     start of a new word.
    for (fst::ArcIterator<Fst> aiter(*fst, s); !aiter.Done(); aiter.Next()) {
      const Arc& arc = aiter.Value();
      if (separators.count(arc.olabel)) {
        new_arcs.push_back(
            Arc(0, 0, fst::Times(arc.weight, bw[arc.nextstate]), sFinal));
        new_arcs_from_init.push_back(
            Arc(0, 0, fst::Times(arc.weight, fw[s]), arc.nextstate));

      } else {
        new_arcs.push_back(arc);
      }
    }
    // Delete all arcs from state s
    fst->DeleteArcs(s);
    // Add new arcs from state s
    for (const Arc& arc : new_arcs) {
      fst->AddArc(s, arc);
    }
  }
  // Add missing arcs from the initial state
  for (const Arc& arc: new_arcs_from_init) {
    fst->AddArc(fst->Start(), arc);
  }
  // Final cost = -total_cost, so that paths are normalized in -log [0, 1]
  fst->SetFinal(sFinal, Weight(-total_cost));
  // Remove epsilon symbols O(V^2 + V * E)
  fst::RmEpsilon(fst);
  // Remove unnecessary states/arcs
  fst::Connect(fst);
  // Empty strings are not allowed
  if (fst->Final(fst->Start()) != Weight::Zero()) {
    fst->SetFinal(fst->Start(), Weight::Zero());
  }
  // Push weights to the toward the initial state. This speeds up n-best list
  // retrieval.
  fst::Push<Arc>(fst, fst::REWEIGHT_TO_INITIAL, fst::kPushWeights);
}

int main(int argc, char** argv) {
  try {
    using namespace kaldi;

    const char* usage =
        "Build a word index from character lattices.\n"
        "\n"
        "The scores of each word are lower bounds of the probability that the "
        "word is present in the transcription, but typically they are very "
        "close to the exact probability.\n"
        "\n"
        "Usage: kaldi-lattice-word-index [options] separator-symbols lat-rspecifier\n"
        " e.g.: kaldi-lattice-word-index \"1 2\" ark:lats.ark\n"
        " e.g.: kaldi-lattice-word-index --nbest=10000 \"1 2\" ark:lats.ark\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    int nbest = 100;
    bool determinize = true;
    bool use_log = false;

    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Register("determinize", &determinize,
                "Determinize fst before building the index.");
    po.Register("nbest", &nbest,
                "Extract this number of paths from the words fst. If the "
                "word fst is deterministic (--determinize=true), then this "
                "is the number of words in the index. Otherwise, the actual "
                "number of words in the index can be smaller.");
    po.Register("use-log", &use_log,
                "If true, perform Forward/Backward using the log semiring "
                "(log-sum-exp), otherwise use the tropical semiring (max).");
    po.Read(argc, argv);

    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    // Parse separator symbols from arguments
    std::unordered_set<fst::StdArc::Label> separator_symbols;
    {
      std::istringstream separator_symbols_iss(po.GetArg(1));
      fst::StdArc::Label tmp;
      while (separator_symbols_iss >> tmp) {
        if (tmp == 0) {
          KALDI_ERR << "Epsilon (0) cannot be a separator symbol!";
        }
        separator_symbols.insert(tmp);
      }
    }

    // Scaling scores
    std::vector<std::vector<double> > scale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = graph_scale;
    scale[1][1] = acoustic_scale;

    const std::string lattice_in_str = po.GetArg(2);
    const bool lattice_is_table =
        (ClassifyRspecifier(lattice_in_str, NULL, NULL) != kNoRspecifier);

    fst::VectorFst<fst::LogArc> log_fst;
    fst::VectorFst<fst::LogArc> pushed_fst;
    fst::VectorFst<fst::StdArc> std_fst;
    fst::VectorFst<fst::StdArc> nbest_fst;

    if (lattice_is_table) {
      SequentialCompactLatticeReader lattice_reader(lattice_in_str);
      for (; !lattice_reader.Done(); lattice_reader.Next()) {
        const std::string lattice_key = lattice_reader.Key();

        CompactLattice lat = lattice_reader.Value();
        lattice_reader.FreeCurrent();
        // Acoustic scale
        if (acoustic_scale != 1.0 || graph_scale != 1.0)
          fst::ScaleLattice(scale, &lat);
        // Word insertion penalty
        if (insertion_penalty != 0.0)
          AddInsPenToLattice(insertion_penalty, &lat);
        // Lattice prunning
        if (beam != std::numeric_limits<BaseFloat>::infinity())
          PruneLattice(beam, &lat);
        // Make sure that lattice complies with all asumptions
        const uint64_t properties =
            lat.Properties(fst::kAcceptor | fst::kAcyclic, true);
        if ((properties & fst::kAcceptor) != fst::kAcceptor) {
          KALDI_ERR << "Lattice " << lattice_key << " is not an acceptor";
        }
        if ((properties & fst::kAcyclic) != fst::kAcyclic) {
          KALDI_ERR << "Lattice " << lattice_key << " is not acyclic";
        }

        if (use_log) {
          // Convert to CompactLattice to VectorFst<LogArc>
          fst::ConvertLattice(lat, &log_fst);
          lat.DeleteStates();
          // Create fst from lattice where each path corresponds to a full
          // WORD in the original lattice.
          WordsFst(&log_fst, separator_symbols);
          // Convert LogArc to StdArc
          fst::ArcMap(log_fst, &std_fst,
                      fst::WeightConvertMapper<fst::LogArc, fst::StdArc>());
          log_fst.DeleteStates();

        } else {
          // Convert to CompactLattice to VectorFst<StdArc>
          fst::ConvertLattice(lat, &std_fst);
          lat.DeleteStates();
          // Create fst from lattice where each path corresponds to a full
          // WORD in the original lattice.
          WordsFst(&std_fst, separator_symbols);
        }

        // Find N-best paths (words)
        if (determinize) {
          fst::ShortestPath(fst::DeterminizeFst<fst::StdArc>(std_fst),
                            &nbest_fst, nbest);
        } else {
          fst::ShortestPath(std_fst, &nbest_fst, nbest);
        }
        std_fst.DeleteStates();

        // Print n-bests and their lower bound score...
        std::vector< fst::VectorFst<fst::StdArc> > nbest_fsts;
        fst::ConvertNbestToVector(nbest_fst, &nbest_fsts);
        for (const fst::VectorFst<fst::StdArc>& fst : nbest_fsts) {
          std::vector<int32> nbest_symbs;
          fst::StdArc::Weight nbest_cost;
          fst::GetLinearSymbolSequence<fst::StdArc, int32>(
              fst, NULL, &nbest_symbs, &nbest_cost);
          std::cout << lattice_key << " " << -nbest_cost.Value();
          for (const int32& s : nbest_symbs)
            std::cout << " " << s;
          std::cout << std::endl;
        }
      }
    } else {
      KALDI_ERR << "NOT IMPLEMENTED!";
    }

    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}
