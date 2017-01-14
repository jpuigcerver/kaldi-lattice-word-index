// MIT License
//
// Copyright (c) 2016 Joan Puigcerver <joapuipe@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

#include <fst/string-weight.h>

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

// Convert a CompactLattice to a FST where the cost of each arc is the
// total cost of the CompactLatticeArc (acoustic + graph costs), the
// input label is the word/char in the CompactLatticeArc, and the output label
// is an integer that maps to the initial and end frames segmentation.
// The mapping table is written in segm_to_label.
template <typename Arc>
int32 CompactLatticeToSegmFst(
    const CompactLattice& clat, fst::MutableFst<Arc>* fst,
    std::map<std::tuple<int32,int32>, int32>* segm_to_label) {
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  fst->DeleteStates();
  segm_to_label->clear();

  // Compute the times for each state in the lattice, so that we can get
  // the segmentation of each symbol.
  std::vector<int32> times;
  const int32 total_frames = CompactLatticeStateTimes(clat, &times);

  // Add states to the output fst.
  for (StateId s = 0; s < clat.NumStates(); ++s) {
    fst->SetFinal(fst->AddState(),
                  clat.Final(s).Weight().Value1() +
                  clat.Final(s).Weight().Value2());
  }
  fst->SetStart(clat.Start());

  // Add arcs to the output fst.
  for (fst::StateIterator<CompactLattice> siter(clat); !siter.Done();
       siter.Next()) {
    const StateId s = siter.Value();
    for (fst::ArcIterator<CompactLattice> aiter(clat, s); !aiter.Done();
         aiter.Next()) {
      const CompactLatticeArc& arc = aiter.Value();
      const std::tuple<int32, int32> segm =
          std::make_tuple(times[s], times[arc.nextstate]);
      const Label new_olabel = segm_to_label->insert(
          std::make_pair(segm, segm_to_label->size() + 1)).first->second;
      const double new_weight =
          arc.weight.Weight().Value1() + arc.weight.Weight().Value2();
      fst->AddArc(s, Arc(arc.ilabel, new_olabel, new_weight, arc.nextstate));
    }
  }

  return total_frames;
}

}  // namespace kaldi

// Given a fst that represents sequences of characters in its paths, creates
// a fst that accepts words (determined by sequences of characters in between
// special labels, separators) which where part of the original fst.
template <typename Arc>
void WordsFst(fst::MutableFst<Arc>* fst,
              const std::unordered_set<int32>& separators) {
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
      if (separators.count(arc.ilabel)) {
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

template <typename Arc>
void PrintIndex(
    const std::string& key, const fst::Fst<Arc>& nbest_fst,
    const std::vector<std::tuple<int32, int32>>& label_to_segm,
    const fst::SymbolTable* symbols_table,
    const bool word_segmentation) {

  // Print n-bests and their scores.
  std::vector< fst::VectorFst<Arc> > nbest_fsts;
  fst::ConvertNbestToVector(nbest_fst, &nbest_fsts);
  for (const fst::VectorFst<Arc>& fst : nbest_fsts) {
    std::vector<int32> nbest_isymbs;
    std::vector<int32> nbest_osymbs;
    typename Arc::Weight nbest_cost;
    fst::GetLinearSymbolSequence<Arc, int32>(
        fst, &nbest_isymbs, &nbest_osymbs, &nbest_cost);
    // Print lattice key and word probability.
    std::cout << key << " " << -nbest_cost.Value();
    // Print char symbols.
    for (const int32& s : nbest_isymbs) {
      if (symbols_table) {
        std::cout << " " << symbols_table->Find(s);
      } else {
        std::cout << " " << s;
      }
    }
    if (word_segmentation) {
      // Print word segmentation: Only first and last frame of the word.
      std::cout << " " << nbest_osymbs.front() - 1
                << " " << nbest_osymbs.back() - 1
                << std::endl;
    } else {
      // Print char segmentation.
      for (const int32& s : nbest_osymbs) {
        // Start frame of each character.
        std::cout << " " << std::get<0>(label_to_segm[s]);
      }
      // Last frame of the last character.
      std::cout << " " << std::get<1>(label_to_segm[nbest_osymbs.back()])
                << std::endl;
    }
  }
}

// Remove the output label from all arcs that are not output arcs from the
// start state and are not inputs to any final state. That is, arcs in between
// of a path.
template <typename Arc>
void RemoveCharSegmFromWordFst(
    fst::MutableFst<Arc>* fst,
    const std::vector<std::tuple<int32,int32>>& label_to_segm) {
  using namespace fst;
  typedef MutableFst<Arc> Fst;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Label Label;
  typedef typename Arc::Weight Weight;

  for (StateIterator<Fst> siter(*fst); !siter.Done(); siter.Next()) {
    const StateId s = siter.Value();
    for (MutableArcIterator<Fst> aiter(fst, s); !aiter.Done(); aiter.Next()) {
      Arc arc = aiter.Value();
      KALDI_ASSERT(arc.olabel > 0);
      if (s == fst->Start()) {
        arc.olabel = std::get<0>(label_to_segm[arc.olabel]) + 1;
      } else if (fst->Final(arc.nextstate) != Weight::Zero()) {
        arc.olabel = std::get<1>(label_to_segm[arc.olabel]) + 1;
      } else {
        arc.olabel = 0;
      }
      aiter.SetValue(arc);
    }
  }
}

int main(int argc, char** argv) {
  try {
    using namespace kaldi;

    const char* usage =
        "Build a word index from character lattices.\n"
        "\n"
        "Words are any sequence of characters in between any of the separator "
        "symbols. The program will output the n-best character segmentations "
        "of words, with their score. More precisely:\n"
        "Let R_c = 1 denote y \\in (\\Sigma^* @)? c_{1:n} (@ \\Sigma^*)? and "
        "R_c = 0 otherwise, where @ is any of the separator symbols.\n"
        "Then R denotes whether the transcription (y) of each utterance "
        "contains the word formed by characters c_1,...,c_n.\n"
        "Let s_{1:n} be a segmentation of each character in c_{1:n}, then "
        "the program computes:\n"
        "\n"
        "1. If --only-best-segmentation=false (the default) then:\n"
        "n-best_{c_{1:n},s_{1:n}} P(R_c = 1| x, s_{1:n})\n"
        "\n"
        "2. If --only-best-segmentation=true then:\n"
        "n-best_{c_{1:n},s_{1:n}} max_{s_{1:n}} P(R_c = 1 | x, s_{1:n})\n"
        "\n"
        "This gives a lower bound to P(R_c = 1 | x), but it is usually quite "
        "close.\n"
        "\n"
        "Usage: kaldi-lattice-word-index [options] separator-symbols "
        "lat-rspecifier\n"
        " e.g.: kaldi-lattice-word-index \"1 2\" ark:lats.ark\n"
        " e.g.: kaldi-lattice-word-index --nbest=10000 \"1 2\" ark:lats.ark\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    BaseFloat delta = fst::kDelta;
    int nbest = 100;
    int32 max_mem = 536870912; // 512MB
    bool only_best_segmentation = false;
    bool word_segmentation = false;
    std::string syms_table_filename = "";

    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Register("delta", &delta, "Tolerance used in determinization. "
                "The smaller the better.");
    po.Register("nbest", &nbest, "Extract this number of n-best hypothesis.");
    po.Register("only-best-segmentation", &only_best_segmentation,
                "If true, output the best character segmentation for each "
                "word.");
    po.Register("word-segmentation", &word_segmentation,
                "If true, output index with the whole word-level segmentation "
                "instead of the character-level segmentation.");
    po.Register("max-mem", &max_mem,
                "Maximum approximate memory usage in determinization (real "
                "usage might be many times this).");
    po.Register("symbols-table", &syms_table_filename,
                "Use this symbols table to map from labels to characters.");
    po.Read(argc, argv);

    if (po.NumArgs() != 2) {
      po.PrintUsage();
      exit(1);
    }

    // Parse separator symbols from arguments
    std::unordered_set<int32> separator_symbols;
    {
      std::istringstream separator_symbols_iss(po.GetArg(1));
      int32 tmp;
      while (separator_symbols_iss >> tmp) {
        if (tmp == 0) {
          KALDI_ERR << "Epsilon (0) cannot be a separator symbol!";
        }
        separator_symbols.insert(tmp);
      }
    }

    // Read symbols table.
    fst::SymbolTable *syms_table = NULL;
    if (syms_table_filename != "") {
      if (!(syms_table = fst::SymbolTable::ReadText(syms_table_filename))) {
        KALDI_ERR << "Could not read symbol table from file "
                  << syms_table_filename;
      }
    }

    // Scaling scores
    std::vector<std::vector<double> > scale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = graph_scale;
    scale[1][1] = acoustic_scale;

    const std::string lattice_in_str = po.GetArg(2);

    fst::VectorFst<fst::LogArc> log_fst;
    std::map<std::tuple<int32, int32>, fst::LogArc::Label> segm_to_label;

    SequentialCompactLatticeReader lattice_reader(lattice_in_str);
    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();

      CompactLattice lat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      // TopSort compact lattice
      TopSortCompactLatticeIfNeeded(&lat);
      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        fst::ScaleLattice(scale, &lat);
      // Word insertion penalty
      if (insertion_penalty != 0.0)
        AddInsPenToLattice(insertion_penalty, &lat);
      // Lattice prunning
      if (beam != std::numeric_limits<BaseFloat>::infinity())
        PruneLattice(beam, &lat);
      // Ensure that the lattice does not have epsilon arcs.
      fst::RmEpsilon(&lat);

      // Convert CompactLattice to LogFst with segmentation info as the
      // output label, and the words/chars as the input label.
      CompactLatticeToSegmFst(lat, &log_fst, &segm_to_label);

      // Create fst from lattice where each path corresponds to a full
      // WORD (and its segmentation) in the original lattice.
      WordsFst(&log_fst, separator_symbols);

      // Reverse the mapping from frame tuples (tuple<int32,int32>) to
      // label (int32).
      std::vector< std::tuple<int32,int32> > label_to_segm;
      label_to_segm.resize(segm_to_label.size() + 1);
      label_to_segm.push_back(std::make_tuple(-1, -1)); // label = 0 not used
      for (const std::pair<std::tuple<int32,int32>, int32>& p : segm_to_label) {
        label_to_segm[p.second] = p.first;
      }

      // We want word-segmentation, instead of the character-level segmentation.
      // Thus, we need to sum all the hypotheses with the same word-level
      // segmentation, even if the character-level segmentation is different.
      // In order to do so, we take advantage of the fact that every path from
      // The initial state to the final state is a word, thus we will remove
      // (set to epsilon) the output label of all arcs except those outgoing from
      // the initial state or entering a final state.
      if (word_segmentation) {
        RemoveCharSegmFromWordFst(&log_fst, label_to_segm);
      }

      // We need to sum up all scores for the same word segmentation.
      // That means determinization in the log-semiring. However, since
      // the FST is non-functional we need to encode the (ilabel, olabel)
      // pairs into a new label and make the original FST a weighted
      // automaton.
      // After determinization, the word and segmentation information are
      // restored again and the weights are mapped to the tropical-semiring
      // in order to get (n-best paths).
      // Note: Determinization, decoding and weight mapping are done on-demand.
      fst::EncodeMapper<fst::LogArc> encoder(fst::kEncodeLabels, fst::ENCODE);
      fst::Encode(&log_fst, &encoder);


      fst::DeterminizeFstOptions<fst::LogArc> det_opts(
          fst::CacheOptions(true, max_mem), delta);
      typedef fst::WeightConvertMapper<fst::LogArc, fst::StdArc> WeightMapper;
      fst::ArcMapFst<fst::LogArc, fst::StdArc, WeightMapper> std_fst(
          fst::DecodeFst<fst::LogArc>(
              fst::DeterminizeFst<fst::LogArc>(log_fst, det_opts), encoder),
          WeightMapper());



      fst::VectorFst<fst::StdArc> nbest_fst;
      if (only_best_segmentation) {
        // We need to determinize again to keep only the best segmentation
        // within the lattice for each word. To do that, we need to
        // determinize again in the tropical-semiring.
        // Also, since the fst is non-functional (for each input sequence
        // (word), we may have multiple output sequences (segmentations).
        fst::DeterminizeFstOptions<fst::StdArc> det_opts2(
            fst::CacheOptions(true, max_mem), delta, 0,
            /* Disambiguate output: */ true);
        fst::ShortestPath(fst::DeterminizeFst<fst::StdArc>(std_fst, det_opts2),
                          &nbest_fst, nbest);
      } else {
        fst::ShortestPath(std_fst, &nbest_fst, nbest);
      }

      // Print index!
      PrintIndex(lattice_key, nbest_fst, label_to_segm, syms_table,
                 word_segmentation);
    }

    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}
